 using namespace std;

#include "BinDef.h"
#include "NewBinDef.h"

#include "iostream"
#include "fstream"

void ReadTheoryFiles(TString FileName, TH1D *h, double scalefactor = 1.0){
    ifstream infile;
    infile.open(FileName.Data());
    double ptmin, ptmax, ptcenter, pp, pperr, AA, AAerr;
    int i = 0;
    string line;
    while (std::getline(infile, line)){
        if (i > 2){ // skip the first 3 lines - They are header information
            std::istringstream iss(line);
            iss >> ptmin >> ptmax >> ptcenter >> pp >> pperr >> AA >> AAerr;
            cout << ptmin << "\t" << ptmax << "\t" << ptcenter << "\t" << pp << "\t" << pperr << "\t" << AA << "\t" << AAerr << endl;
            h->SetBinContent(i-2, AA);
            h->SetBinError(i-2, AAerr);
            // cout << AA << "\t" << AAerr << endl;
        }
        i++;
    }
    h->Scale(scalefactor);    
}

void ReadRCPTheoryFiles(TString FileName, TH1D *h, double scalefactor = 1.0){
    ifstream infile;
    infile.open(FileName.Data());
    double ptmin, ptmax, ptcenter, pp, pperr, AA, AAerr;
    int i = 0;
    string line;
    while (std::getline(infile, line)){
        if (i > 2){ // skip the first 3 lines - They are header information
            std::istringstream iss(line);
            iss >> ptmin >> ptmax >> ptcenter >> AA >> AAerr;
            cout << ptmin << "\t" << ptmax << "\t" << ptcenter << "\t" << AA << "\t" << AAerr << endl;
            h->SetBinContent(i-2, AA);
            h->SetBinError(i-2, AAerr);
            // cout << AA << "\t" << AAerr << endl;
        }
        i++;
    }
    h->Scale(scalefactor);    
}

void SetBinError(TH1 *h1, TH1 *h2){
    assert(h1->GetNbinsX() == h2->GetNbinsX());
    for (int i = 1; i <= h1->GetNbinsX(); i++){
        h1->SetBinError(i, h2->GetBinError(i));
    }
}

// TGraphAsymmErrors *GraphPlots(TH1 *Main, vector<TH1 *>low, vector<TH1 *>high){
//     assert(low.size() == high.size());
//     TH1D *TotalLow = (TH1D *)low[0]->Clone();
//     TotalLow->Reset();
//     for (int j = 1; j <= low[i]->GetNbinsX(); j++){
//         double totalerr = 0;
//         for (int i = 0; i < low.size(); i++){     
//             totalerr += pow(low[i]->GetBinError(j), 2);
//         }
//         totalerr = sqrt(totalerr);
//         TotalLow->SetBinContent(j, totalerr);
//     }

//     TH1D *TotalHigh = (TH1D *)high[0]->Clone();
//     TotalHigh->Reset();
//     for (int j = 1; j <= high[i]->GetNbinsX(); j++){
//         double totalerr = 0;
//         for (int i = 0; i < high.size(); i++){     
//             totalerr += pow(high[i]->GetBinError(j), 2);
//         }
//         totalerr = sqrt(totalerr);
//         TotalHigh->SetBinContent(j, totalerr);
//     }

//     TGraphAsymmErrors *g = new TGraphAsymmErrors(Main);
//     for (int j = 1; j <= low[i]->GetNbinsX(); j++){

//         g->SetPointEYlow(j-1, 0.01*TotalLow->GetBinError(j)*Main->GetBinContent(j));
//         g->SetPointEYhigh(j-1, 0.01*TotalHigh->GetBinError(j)*Main->GetBinContent(j));
//     }
//     return g;
// }


vector<TString> BaseDirName;
vector<TString> ExtraTag;
vector<int> nPtBins;
vector<int> nZBins;
vector<int> nIter;

void AddASet(TString BaseDir, TString Extra, int PtBins, int ZBins, int iter){
    ::BaseDirName.push_back(BaseDir);
    ::ExtraTag.push_back(Extra);
    ::nPtBins.push_back(PtBins);
    ::nZBins.push_back(ZBins);
    ::nIter.push_back(iter);
}

int color[10] = {kRed, kGreen-2, kBlue, kBlack, kCyan, kMagenta, kOrange, kOrange+2, 5, 5};
int markerstyle[10] = { kOpenCircle, kOpenSquare, kOpenTriangleUp, kOpenTriangleDown, kFullTriangleUp, kFullTriangleDown, kFullSquare, kFullCircle, 29, 30};

double taa[3] = {941.23714, 391.35550, 56.62475};
double nevents[3] = {1.0318440e+08, 3.2123506e+08, 4.6679240e+08};

TString JetPtXaxisName = "p_{T, Jet}^{D^{0}} [GeV/#it{c}]";
TString JetZXaxisName = "z_{Jet}^{D^{0}}";
TString JetdRXaxisName = "#DeltaR";

TString JetPtYaxisName = "#frac{1}{N_{Evt}} #frac{d^{2}N}{p_{T}dp_{T, Jet}d#eta} [GeV/#it{c}]^{-2}";
TString JetZYaxisName = "#frac{1}{N_{Evt}} #frac{d^{2}N}{zdz_{Jet}d#eta}";
TString JetdRYaxisName = "#frac{1}{N_{Jet}} #frac{#DeltaN_{Jet}}{#DeltaR}";
TString SystematicUncYaxisName = "(Variation - Nominal)/Nominal";

TString RCPYaxisName = "R_{CP}";
TString RMPYaxisName = "R_{MP}";

void CompareSpectraAndRCPBetweenMethods(){
    TH1::SetDefaultSumw2();
  	TH2::SetDefaultSumw2();
  	gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    TString PlotDir = "Plots";
    bool use1Dhistogram = false;
    if (use1Dhistogram) PlotDir += "_1D";

    AddASet("Aug14_FONLL", "Data", 12, 7, 4);
    // Prior Variation Systematic
    AddASet("Aug14_FONLL", "Data_WeighedByData", 12, 7, 4);
    AddASet("Aug14_PYTHIA", "Data", 12, 7, 4);
    // Iteration Parameter Systematic
    AddASet("Aug14_FONLL", "Data", 12, 7, 3);
    AddASet("Aug14_FONLL", "Data", 12, 7, 5);
    // Tracking Efficiency Systematic
    AddASet("Aug14_FONLL_TrackEff", "Data", 12, 7, 4);
    // Signal Extraction Systematic
    AddASet("Aug14_FONLL_D0YieldLower", "Data", 12, 7, 4);
    AddASet("Aug14_FONLL_D0YieldUpper", "Data", 12, 7, 4);
    // D0 Reconstruction Efficiency Systematic
    AddASet("Aug14_FONLL_D0YieldEffLower", "Data", 12, 7, 4);
    AddASet("Aug14_FONLL_D0YieldEffUpper", "Data", 12, 7, 4);


    if (BaseDirName.size() != ExtraTag.size()){
        cout << "BaseDirName and ExtraTag must be the same size!" << endl;
        return;
    }
    // const int nD0PtBins = 5;
    const int nCent = 3;
    const int nSize = BaseDirName.size();
    const int nType = 3; // Types of Model Variation

    TH1D *AvgJetpT[5][3];
    TH1D *AvgJetZ[5][3];
    TH1D *AvgJetZPerJet[5][3];
    TH1D *AvgJetdR[5][3];

    TH1D *AvgRCP_JetpT[5][3];
    TH1D *AvgRCP_JetZ[5][3];
    TH1D *AvgRCP_JetZPerJet[5][3];
    TH1D *AvgRCP_JetdR[5][3];

    TH1D *JetpT[5][100][3];
    TH1D *JetZ[5][100][3];
    TH1D *JetZPerJet[5][100][3];
    TH1D *JetdR[5][100][3];
    TH1D *RCP_Pt[5][100][3];
    TH1D *RCP_Z[5][100][3];
    TH1D *RCP_ZPerJet[5][100][3];
    TH1D *RCP_dR[5][100][3];

    // Systematic Uncertainty From Iteration From Unfolding, Prior From PYTHIA
    
    vector<TString> SystematicName = {"Iter=3", "Iter=5", "Data-Weighed", "PYTHIA-Fit", "Tracking-Efficiency", "D0-Yield-Lower", "D0-Yield-Upper", "D0-Efficiency-Lower", "D0-Efficiency-Upper"};
    vector<int> SystematicHist = {3, 4, 1, 2, 5, 6, 7, 8, 9};

    const int nSystematic = SystematicHist.size();

    TH1D *SystematicUncertainty_Pt[5][nSystematic][3];
    TH1D *SystematicUncertainty_Z[5][nSystematic][3];
    TH1D *SystematicUncertainty_ZPerJet[5][nSystematic][3];
    TH1D *SystematicUncertainty_dR[5][nSystematic][3];

    TH1D *SystematicUncertainty_RCP_Pt[5][nSystematic][2];
    TH1D *SystematicUncertainty_RCP_Z[5][nSystematic][2];
    TH1D *SystematicUncertainty_RCP_ZPerJet[5][nSystematic][2];
    TH1D *SystematicUncertainty_RCP_dR[5][nSystematic][2];

    cout << "Assigning Histograms" << endl;

    for (int pT = 1; pT <= 5; pT++){
        for (int i = 0; i < nSize; i++){   
            TFile *f = new TFile(Form("%s_%iGeV_%s%i/Output_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", BaseDirName[i].Data(), pT, ExtraTag[i].Data(), nIter[i], 0, 0, nPtBins[i], nZBins[i]));
            cout << "FileName is " << f->GetName() << endl;
            f->cd();
            for (int cent = 0; cent < 3; cent++){
                if (!use1Dhistogram) JetpT[pT-1][i][cent] = (TH1D *)gDirectory->Get(Form("Unfolded Wide p_{T} Cent = %i SI = %i Iter = %i", cent, 0, nIter[i]));
                else JetpT[pT-1][i][cent] = (TH1D *)gDirectory->Get(Form("Unfolded1D Cent = %i SI = %i Iter = %i", cent, 0, nIter[i]));
                cout << JetpT[pT-1][i][cent]->GetName() << endl;
                SetName(JetpT[pT-1][i][cent], Form("%s %s D0pT = %i GeV", JetpT[pT-1][i][cent]->GetName(), ExtraTag[i].Data(), pT));
                SetAxisTitles(JetpT[pT-1][i][cent], JetPtXaxisName.Data(), JetPtYaxisName.Data());
                JetZ[pT-1][i][cent] = (TH1D *)gDirectory->Get(Form("Unfolded Wide Z Cent = %i SI = %i Iter = %i", cent, 0, nIter[i]));
                // cout << JetZ[i][cent]->GetName() << endl;
                SetName(JetZ[pT-1][i][cent], Form("%s %s D0pT = %i GeV", JetZ[pT-1][i][cent]->GetName(), ExtraTag[i].Data(), pT));
                SetAxisTitles(JetZ[pT-1][i][cent], JetZXaxisName.Data(), JetZYaxisName.Data());
                JetZPerJet[pT-1][i][cent] = (TH1D *)gDirectory->Get(Form("Unfolded Wide Z Cent = %i SI = %i Iter = %i", cent, 0, nIter[i]));
                // cout << JetZ[i][cent]->GetName() << endl;
                SetName(JetZPerJet[pT-1][i][cent], Form("%s Per Jet %s D0pT = %i GeV", JetZ[pT-1][i][cent]->GetName(), ExtraTag[i].Data(), pT));
                SetAxisTitles(JetZPerJet[pT-1][i][cent], JetZXaxisName.Data(), JetZYaxisName.Data());
                JetdR[pT-1][i][cent] = (TH1D *)gDirectory->Get(Form("Unfolded Wide dR Cent = %i SI = %i Iter = %i", cent, 0, nIter[i]));
                // cout << JetdR[i][cent]->GetName() << endl;
                SetName(JetdR[pT-1][i][cent], Form("%s %s D0pT = %i GeV", JetdR[pT-1][i][cent]->GetName(), ExtraTag[i].Data(), pT));
                SetAxisTitles(JetdR[pT-1][i][cent], JetdRXaxisName.Data(), JetdRYaxisName.Data());
                // cout << "================ Jet dR R = 0.2 to R = 0.4 counts ================" << JetdR[pT-1][i][cent]->GetBinContent(4)/JetdR[pT-1][i][cent]->Integral() << endl;
                JetpT[pT-1][i][cent]->Scale(1./nevents[cent]);
                JetZ[pT-1][i][cent]->Scale(1./nevents[cent]);
                ProcessSpectra(JetpT[pT-1][i][cent]);
                ProcessSpectra(JetZ[pT-1][i][cent]);
                JetZPerJet[pT-1][i][cent]->Scale(1./JetZPerJet[pT-1][i][cent]->Integral());
                DivideByBinWidth(JetZPerJet[pT-1][i][cent]);
                JetdR[pT-1][i][cent]->Scale(1./JetdR[pT-1][i][cent]->Integral());
                DivideByBinWidth(JetdR[pT-1][i][cent]);
                // cout << "Done" << endl;
                // Color scheme: Different set of iterations have different colors. The markers are consistent among FONLL, FONLL_DataWeighed, PYTHIA
                const int type = i % nType; // This is the prior type. Useful only for plotting.
                SetColor(JetpT[pT-1][i][cent], color[i], markerstyle[i], 1);
                SetColor(JetZ[pT-1][i][cent], color[i], markerstyle[i], 1);
                SetColor(JetZPerJet[pT-1][i][cent], color[i], markerstyle[i], 1);
                SetColor(JetdR[pT-1][i][cent], color[i], markerstyle[i], 1);

                JetpT[pT-1][i][cent]->SetDirectory(0);
                JetZ[pT-1][i][cent]->SetDirectory(0);
                JetZPerJet[pT-1][i][cent]->SetDirectory(0);
                JetdR[pT-1][i][cent]->SetDirectory(0);

                cout << JetpT[pT-1][i][cent]->GetNbinsX() << "\t" << JetZ[pT-1][i][cent]->GetNbinsX() << "\t" << JetdR[pT-1][i][cent]->GetNbinsX() << endl;

                RCP_Pt[pT-1][i][cent] = (TH1D *)JetpT[pT-1][i][cent]->Clone();
                SetName(RCP_Pt[pT-1][i][cent], Form("%s RCP", JetpT[pT-1][i][cent]->GetName()));
                RCP_Z[pT-1][i][cent] = (TH1D *)JetZ[pT-1][i][cent]->Clone();
                SetName(RCP_Z[pT-1][i][cent], Form("%s RCP", JetZ[pT-1][i][cent]->GetName()));
                RCP_ZPerJet[pT-1][i][cent] = (TH1D *)JetZPerJet[pT-1][i][cent]->Clone();
                SetName(RCP_ZPerJet[pT-1][i][cent], Form("%s RCP", JetZPerJet[pT-1][i][cent]->GetName()));
                RCP_dR[pT-1][i][cent] = (TH1D *)JetdR[pT-1][i][cent]->Clone();
                SetName(RCP_dR[pT-1][i][cent], Form("%s RCP", JetdR[pT-1][i][cent]->GetName()));

                RCP_Pt[pT-1][i][cent]->SetDirectory(0);
                RCP_Z[pT-1][i][cent]->SetDirectory(0);
                RCP_ZPerJet[pT-1][i][cent]->SetDirectory(0);
                RCP_dR[pT-1][i][cent]->SetDirectory(0);
            }
            f->Close();
            

            RCP_Pt[pT-1][i][0]->Divide(RCP_Pt[pT-1][i][2]);
            RCP_Pt[pT-1][i][0]->Scale(taa[2]/taa[0]);
            RCP_Pt[pT-1][i][1]->Divide(RCP_Pt[pT-1][i][2]);
            RCP_Pt[pT-1][i][1]->Scale(taa[2]/taa[1]);
            RCP_Z[pT-1][i][0]->Divide(RCP_Z[pT-1][i][2]);
            RCP_Z[pT-1][i][0]->Scale(taa[2]/taa[0]);
            RCP_Z[pT-1][i][1]->Divide(RCP_Z[pT-1][i][2]);
            RCP_Z[pT-1][i][1]->Scale(taa[2]/taa[1]);
            RCP_ZPerJet[pT-1][i][0]->Divide(RCP_ZPerJet[pT-1][i][2]);
            RCP_ZPerJet[pT-1][i][1]->Divide(RCP_ZPerJet[pT-1][i][2]);
            RCP_dR[pT-1][i][0]->Divide(RCP_dR[pT-1][i][2]);
            RCP_dR[pT-1][i][1]->Divide(RCP_dR[pT-1][i][2]);     
            
        }
    }

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){   
            AvgJetpT[pT-1][cent] = (TH1D *)JetpT[pT-1][0][cent]->Clone();
            SetName(AvgJetpT[pT-1][cent], Form("%s Avg D0pT = %i GeV", JetpT[pT-1][0][cent]->GetName(), pT));
            AvgJetpT[pT-1][cent]->Add(JetpT[pT-1][1][cent]);
            AvgJetpT[pT-1][cent]->Add(JetpT[pT-1][2][cent]);
            AvgJetpT[pT-1][cent]->Scale(1./3.);

            SetBinError(AvgJetpT[pT-1][cent], JetpT[pT-1][0][cent]);

            AvgJetZ[pT-1][cent] = (TH1D *)JetZ[pT-1][0][cent]->Clone();
            SetName(AvgJetZ[pT-1][cent], Form("%s Avg D0pT = %i GeV", JetZ[pT-1][0][cent]->GetName(), pT));
            AvgJetZ[pT-1][cent]->Add(JetZ[pT-1][1][cent]);
            AvgJetZ[pT-1][cent]->Add(JetZ[pT-1][2][cent]);
            AvgJetZ[pT-1][cent]->Scale(1./3.);

            SetBinError(AvgJetZ[pT-1][cent], JetZ[pT-1][0][cent]);

            AvgJetZPerJet[pT-1][cent] = (TH1D *)JetZPerJet[pT-1][0][cent]->Clone();
            SetName(AvgJetZPerJet[pT-1][cent], Form("%s Avg D0pT = %i GeV", JetZPerJet[pT-1][0][cent]->GetName(), pT));
            AvgJetZPerJet[pT-1][cent]->Add(JetZPerJet[pT-1][1][cent]);
            AvgJetZPerJet[pT-1][cent]->Add(JetZPerJet[pT-1][2][cent]);
            AvgJetZPerJet[pT-1][cent]->Scale(1./3.);

            SetBinError(AvgJetZPerJet[pT-1][cent], JetZPerJet[pT-1][0][cent]);

            AvgJetdR[pT-1][cent] = (TH1D *)JetdR[pT-1][0][cent]->Clone();
            SetName(AvgJetdR[pT-1][cent], Form("%s Avg D0pT = %i GeV", JetdR[pT-1][0][cent]->GetName(), pT));
            AvgJetdR[pT-1][cent]->Add(JetdR[pT-1][1][cent]);
            AvgJetdR[pT-1][cent]->Add(JetdR[pT-1][2][cent]);
            AvgJetdR[pT-1][cent]->Scale(1./3.);

            SetBinError(AvgJetdR[pT-1][cent], JetdR[pT-1][0][cent]);

            AvgRCP_JetpT[pT-1][cent] = (TH1D *)RCP_Pt[pT-1][0][cent]->Clone();
            SetName(AvgRCP_JetpT[pT-1][cent], Form("%s Avg D0pT = %i GeV", RCP_Pt[pT-1][0][cent]->GetName(), pT));
            AvgRCP_JetpT[pT-1][cent]->Add(RCP_Pt[pT-1][1][cent]);
            AvgRCP_JetpT[pT-1][cent]->Add(RCP_Pt[pT-1][2][cent]);
            AvgRCP_JetpT[pT-1][cent]->Scale(1./3.);

            SetBinError(AvgRCP_JetpT[pT-1][cent], RCP_Pt[pT-1][0][cent]);

            AvgRCP_JetZ[pT-1][cent] = (TH1D *)RCP_Z[pT-1][0][cent]->Clone();
            SetName(AvgRCP_JetZ[pT-1][cent], Form("%s Avg D0pT = %i GeV", RCP_Z[pT-1][0][cent]->GetName(), pT));
            AvgRCP_JetZ[pT-1][cent]->Add(RCP_Z[pT-1][1][cent]);
            AvgRCP_JetZ[pT-1][cent]->Add(RCP_Z[pT-1][2][cent]);
            AvgRCP_JetZ[pT-1][cent]->Scale(1./3.);

            SetBinError(AvgRCP_JetZ[pT-1][cent], RCP_Z[pT-1][0][cent]);

            AvgRCP_JetZPerJet[pT-1][cent] = (TH1D *)RCP_ZPerJet[pT-1][0][cent]->Clone();
            SetName(AvgRCP_JetZPerJet[pT-1][cent], Form("%s Avg D0pT = %i GeV", RCP_ZPerJet[pT-1][0][cent]->GetName(), pT));
            AvgRCP_JetZPerJet[pT-1][cent]->Add(RCP_ZPerJet[pT-1][1][cent]);
            AvgRCP_JetZPerJet[pT-1][cent]->Add(RCP_ZPerJet[pT-1][2][cent]);
            AvgRCP_JetZPerJet[pT-1][cent]->Scale(1./3.);

            SetBinError(AvgRCP_JetZPerJet[pT-1][cent], RCP_ZPerJet[pT-1][0][cent]);

            AvgRCP_JetdR[pT-1][cent] = (TH1D *)RCP_dR[pT-1][0][cent]->Clone();
            SetName(AvgRCP_JetdR[pT-1][cent], Form("%s Avg D0pT = %i GeV", RCP_dR[pT-1][0][cent]->GetName(), pT));
            AvgRCP_JetdR[pT-1][cent]->Add(RCP_dR[pT-1][1][cent]);
            AvgRCP_JetdR[pT-1][cent]->Add(RCP_dR[pT-1][2][cent]);
            AvgRCP_JetdR[pT-1][cent]->Scale(1./3.);

            SetBinError(AvgRCP_JetdR[pT-1][cent], RCP_dR[pT-1][0][cent]);
        }
    }

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            for (int i = 0; i < nSystematic; i++){
                SystematicUncertainty_Pt[pT-1][i][cent] = (TH1D *)JetpT[pT-1][SystematicHist[i]][cent]->Clone();
                SetName(SystematicUncertainty_Pt[pT-1][i][cent], Form("%s Systematic Unc. %s", JetpT[pT-1][SystematicHist[i]][cent]->GetName(), SystematicName[i].Data()));
                cout << SystematicUncertainty_Pt[pT-1][i][cent]->GetName() << "\t" << SystematicUncertainty_Pt[pT-1][i][cent]->GetNbinsX() << endl;
                SystematicUncertainty_Pt[pT-1][i][cent]->Add(JetpT[pT-1][0][cent], -1);
                SystematicUncertainty_Pt[pT-1][i][cent]->Divide(JetpT[pT-1][0][cent]);
                SystematicUncertainty_Pt[pT-1][i][cent]->Scale(100.);
            }

            for (int i = 0; i < nSystematic; i++){
                SystematicUncertainty_Z[pT-1][i][cent] = (TH1D *)JetZ[pT-1][SystematicHist[i]][cent]->Clone();
                SetName(SystematicUncertainty_Z[pT-1][i][cent], Form("%s Systematic Unc. %s", JetZ[pT-1][SystematicHist[i]][cent]->GetName(), SystematicName[i].Data()));
                SystematicUncertainty_Z[pT-1][i][cent]->Add(JetZ[pT-1][0][cent], -1);
                SystematicUncertainty_Z[pT-1][i][cent]->Divide(JetZ[pT-1][0][cent]);
                SystematicUncertainty_Z[pT-1][i][cent]->Scale(100.);
            }

            for (int i = 0; i < nSystematic; i++){
                SystematicUncertainty_ZPerJet[pT-1][i][cent] = (TH1D *)JetZPerJet[pT-1][SystematicHist[i]][cent]->Clone();
                SetName(SystematicUncertainty_ZPerJet[pT-1][i][cent], Form("%s Systematic Unc. %s", JetZPerJet[pT-1][SystematicHist[i]][cent]->GetName(), SystematicName[i].Data()));
                SystematicUncertainty_ZPerJet[pT-1][i][cent]->Add(JetZPerJet[pT-1][0][cent], -1);
                SystematicUncertainty_ZPerJet[pT-1][i][cent]->Divide(JetZPerJet[pT-1][0][cent]);
                SystematicUncertainty_ZPerJet[pT-1][i][cent]->Scale(100.);
            }

            for (int i = 0; i < nSystematic; i++){
                SystematicUncertainty_dR[pT-1][i][cent] = (TH1D *)JetdR[pT-1][SystematicHist[i]][cent]->Clone();
                SetName(SystematicUncertainty_dR[pT-1][i][cent], Form("%s Systematic Unc. %s", JetdR[pT-1][SystematicHist[i]][cent]->GetName(), SystematicName[i].Data()));
                SystematicUncertainty_dR[pT-1][i][cent]->Add(JetdR[pT-1][0][cent], -1);
                SystematicUncertainty_dR[pT-1][i][cent]->Divide(JetdR[pT-1][0][cent]);
                SystematicUncertainty_dR[pT-1][i][cent]->Scale(100.);
            }
        }
    }

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){
            for (int i = 0; i < nSystematic; i++){
                SystematicUncertainty_RCP_Pt[pT-1][i][cent] = (TH1D *)RCP_Pt[pT-1][SystematicHist[i]][cent]->Clone();
                SetName(SystematicUncertainty_RCP_Pt[pT-1][i][cent], Form("%s Systematic Unc. %s", RCP_Pt[pT-1][SystematicHist[i]][cent]->GetName(), SystematicName[i].Data()));
                cout << SystematicUncertainty_RCP_Pt[pT-1][i][cent]->GetName() << "\t" << SystematicUncertainty_RCP_Pt[pT-1][i][cent]->GetNbinsX() << endl;
                SystematicUncertainty_RCP_Pt[pT-1][i][cent]->Add(RCP_Pt[pT-1][0][cent], -1);
                SystematicUncertainty_RCP_Pt[pT-1][i][cent]->Divide(RCP_Pt[pT-1][0][cent]);
                SystematicUncertainty_RCP_Pt[pT-1][i][cent]->Scale(100.);
            }

            for (int i = 0; i < nSystematic; i++){
                SystematicUncertainty_RCP_Z[pT-1][i][cent] = (TH1D *)RCP_Z[pT-1][SystematicHist[i]][cent]->Clone();
                SetName(SystematicUncertainty_RCP_Z[pT-1][i][cent], Form("%s Systematic Unc. %s", RCP_Z[pT-1][SystematicHist[i]][cent]->GetName(), SystematicName[i].Data()));
                SystematicUncertainty_RCP_Z[pT-1][i][cent]->Add(RCP_Z[pT-1][0][cent], -1);
                SystematicUncertainty_RCP_Z[pT-1][i][cent]->Divide(RCP_Z[pT-1][0][cent]);
                SystematicUncertainty_RCP_Z[pT-1][i][cent]->Scale(100.);
            }

            for (int i = 0; i < nSystematic; i++){
                SystematicUncertainty_RCP_ZPerJet[pT-1][i][cent] = (TH1D *)RCP_ZPerJet[pT-1][SystematicHist[i]][cent]->Clone();
                SetName(SystematicUncertainty_RCP_ZPerJet[pT-1][i][cent], Form("%s Systematic Unc. %s", RCP_ZPerJet[pT-1][SystematicHist[i]][cent]->GetName(), SystematicName[i].Data()));
                SystematicUncertainty_RCP_ZPerJet[pT-1][i][cent]->Add(RCP_ZPerJet[pT-1][0][cent], -1);
                SystematicUncertainty_RCP_ZPerJet[pT-1][i][cent]->Divide(RCP_ZPerJet[pT-1][0][cent]);
                SystematicUncertainty_RCP_ZPerJet[pT-1][i][cent]->Scale(100.);
            }

            for (int i = 0; i < nSystematic; i++){
                SystematicUncertainty_RCP_dR[pT-1][i][cent] = (TH1D *)RCP_dR[pT-1][SystematicHist[i]][cent]->Clone();
                SetName(SystematicUncertainty_RCP_dR[pT-1][i][cent], Form("%s Systematic Unc. %s", RCP_dR[pT-1][SystematicHist[i]][cent]->GetName(), SystematicName[i].Data()));
                SystematicUncertainty_RCP_dR[pT-1][i][cent]->Add(RCP_dR[pT-1][0][cent], -1);
                SystematicUncertainty_RCP_dR[pT-1][i][cent]->Divide(RCP_dR[pT-1][0][cent]);
                SystematicUncertainty_RCP_dR[pT-1][i][cent]->Scale(100.);
            }
        }
    }
    
    TCanvas *PtCanvas = new TCanvas("PtCanvas", "PtCanvas", 2100, 700);
    PtCanvas->Divide(5,3);

    TLegend *legend[5];
    
    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            if (cent == 0) legend[pT-1] = new TLegend(0.1,0.1,0.48,0.3);
            PtCanvas->cd(cent*5+pT);
            gPad->SetLogy();
            for (int i = 0; i < nSize; i++){
                JetpT[pT-1][i][cent]->Draw("SAME");
                JetpT[pT-1][i][cent]->GetXaxis()->SetRangeUser(5,20);
                JetpT[pT-1][i][cent]->GetYaxis()->SetRangeUser(1e-10,1e-3);
                TString LegendEntry = BaseDirName[i]+ExtraTag[i]+TString(to_string(nIter[i]));
                if (cent == 0) legend[pT-1]->AddEntry(JetpT[pT-1][i][cent], LegendEntry.Data(), "lp");
            }
            legend[pT-1]->Draw("SAME");
        }
        
    }

    PtCanvas->SaveAs(Form("%s/Pt_Iter_%i.pdf", PlotDir.Data(),4));

    TCanvas *ZCanvas = new TCanvas("ZCanvas", "ZCanvas", 2100, 700);
    ZCanvas->Divide(5,3);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            ZCanvas->cd(cent*5+pT);
            gPad->SetLogy();
            for (int i = 0; i < nSize; i++){
                JetZ[pT-1][i][cent]->Draw("SAME");
                JetZ[pT-1][i][cent]->GetYaxis()->SetRangeUser(1e-8,1e-1);
            }
            legend[pT-1]->Draw("SAME");
        }
        
    }

    ZCanvas->SaveAs(Form("%s/Z_Iter_%i.pdf", PlotDir.Data(),4));

    TCanvas *ZPerJetCanvas = new TCanvas("ZPerJetCanvas", "ZPerJetCanvas", 2100, 700);
    ZPerJetCanvas->Divide(5,3);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            ZPerJetCanvas->cd(cent*5+pT);
            gPad->SetLogy();
            for (int i = 0; i < nSize; i++){
                JetZPerJet[pT-1][i][cent]->Draw("SAME");
                JetZPerJet[pT-1][i][cent]->GetYaxis()->SetRangeUser(1e-4,1e4);
            }
            legend[pT-1]->Draw("SAME");
        }
        
    }

    ZPerJetCanvas->SaveAs(Form("%s/ZPerJet_Iter_%i.pdf", PlotDir.Data(),4));

    TCanvas *dRCanvas = new TCanvas("dRCanvas", "dRCanvas", 2100, 700);
    dRCanvas->Divide(5,3);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            dRCanvas->cd(cent*5+pT);
            gPad->SetLogy();
            for (int i = 0; i < nSize; i++){
                JetdR[pT-1][i][cent]->Draw("SAME");
                JetdR[pT-1][i][cent]->GetXaxis()->SetRangeUser(0, 0.2);
                JetdR[pT-1][i][cent]->GetYaxis()->SetRangeUser(1e-4,1e4);
            }
            legend[pT-1]->Draw("SAME");
        }
        
    }

    dRCanvas->SaveAs(Form("%s/dR_Iter_%i.pdf", PlotDir.Data(),4));


    TCanvas *RCP_PtCanvas = new TCanvas("RCP_PtCanvas", "RCP_PtCanvas", 2100, 700);
    RCP_PtCanvas->Divide(5,2);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){
            RCP_PtCanvas->cd(cent*5+pT);
            for (int i = 0; i < nSize; i++){
                RCP_Pt[pT-1][i][cent]->GetYaxis()->SetRangeUser(0,2);
                if (cent == 0) RCP_Pt[pT-1][i][cent]->GetYaxis()->SetTitle(RCPYaxisName.Data());
                else RCP_Pt[pT-1][i][cent]->GetYaxis()->SetTitle(RMPYaxisName.Data());
                RCP_Pt[pT-1][i][cent]->GetXaxis()->SetRangeUser(5,20);
                RCP_Pt[pT-1][i][cent]->Draw("SAME");
            }
            legend[pT-1]->Draw("SAME");
        }
        
    }

    RCP_PtCanvas->SaveAs(Form("%s/RCP_Pt_Iter_%i.pdf", PlotDir.Data(),4));


    TCanvas *RCP_ZCanvas = new TCanvas("RCP_ZCanvas", "RCP_ZCanvas", 2100, 700);
    RCP_ZCanvas->Divide(5,2);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){
            RCP_ZCanvas->cd(cent*5+pT);
            for (int i = 0; i < nSize; i++){
                cout << "Drawing " << RCP_Z[pT-1][i][cent]->GetName() << endl;
                RCP_Z[pT-1][i][cent]->GetYaxis()->SetRangeUser(0,2);
                if (cent == 0) RCP_Z[pT-1][i][cent]->GetYaxis()->SetTitle(RCPYaxisName.Data());
                else RCP_Z[pT-1][i][cent]->GetYaxis()->SetTitle(RMPYaxisName.Data());
                RCP_Z[pT-1][i][cent]->Draw("SAME");
            }
            legend[pT-1]->Draw("SAME");
        }
        
    }

    RCP_ZCanvas->SaveAs(Form("%s/RCP_Z_Iter_%i.pdf", PlotDir.Data(),4));

    TCanvas *RCP_ZPerJetCanvas = new TCanvas("RCP_ZPerJetCanvas", "RCP_ZPerJetCanvas", 2100, 700);
    RCP_ZPerJetCanvas->Divide(5,2);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){
            RCP_ZPerJetCanvas->cd(cent*5+pT);
            for (int i = 0; i < nSize; i++){
                cout << "Drawing " << RCP_ZPerJet[pT-1][i][cent]->GetName() << endl;
                RCP_ZPerJet[pT-1][i][cent]->GetYaxis()->SetRangeUser(0,2);
                if (cent == 0) RCP_ZPerJet[pT-1][i][cent]->GetYaxis()->SetTitle(RCPYaxisName.Data());
                else RCP_ZPerJet[pT-1][i][cent]->GetYaxis()->SetTitle(RMPYaxisName.Data());
                RCP_ZPerJet[pT-1][i][cent]->Draw("SAME");
            }
            legend[pT-1]->Draw("SAME");
        }
        
    }

    RCP_ZPerJetCanvas->SaveAs(Form("%s/RCP_ZPerJet_Iter_%i.pdf", PlotDir.Data(),4));

    TCanvas *RCP_dRCanvas = new TCanvas("RCP_dRCanvas", "RCP_dRCanvas", 2100, 700);
    RCP_dRCanvas->Divide(5,2);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){
            RCP_dRCanvas->cd(cent*5+pT);
            for (int i = 0; i < nSize; i++){
                cout << "Drawing " << RCP_dR[pT-1][i][cent]->GetName() << endl;
                RCP_dR[pT-1][i][cent]->GetYaxis()->SetRangeUser(0.5,1.5);
                RCP_dR[pT-1][i][cent]->GetXaxis()->SetRangeUser(0,0.2);
                if (cent == 0) RCP_dR[pT-1][i][cent]->GetYaxis()->SetTitle(RCPYaxisName.Data());
                else RCP_dR[pT-1][i][cent]->GetYaxis()->SetTitle(RMPYaxisName.Data());
                RCP_dR[pT-1][i][cent]->Draw("SAME");
            }
            legend[pT-1]->Draw("SAME");
        }
        
    }

    RCP_dRCanvas->SaveAs(Form("%s/RCP_dR_Iter_%i.pdf", PlotDir.Data(),4));

    auto sysunclegend = new TLegend(0.1,0.1,0.48,0.4);

    TCanvas *SystematicUncertainty_PtCanvas = new TCanvas("SystematicUncertainty_PtCanvas", "SystematicUncertainty_PtCanvas", 2100, 700);
    SystematicUncertainty_PtCanvas->Divide(5,3);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            SystematicUncertainty_PtCanvas->cd(cent*5+pT);
            for (int i = 0; i < nSystematic; i++){
                SystematicUncertainty_Pt[pT-1][i][cent]->GetYaxis()->SetTitle(SystematicUncYaxisName.Data());
                SystematicUncertainty_Pt[pT-1][i][cent]->GetYaxis()->SetRangeUser(-100,100);
                SystematicUncertainty_Pt[pT-1][i][cent]->GetXaxis()->SetRangeUser(5,20);
                SystematicUncertainty_Pt[pT-1][i][cent]->Draw("LP HIST SAME");
                if (cent == 0 && pT == 1) sysunclegend->AddEntry(SystematicUncertainty_Pt[pT-1][i][cent], SystematicName[i].Data(), "lp");
            }
            sysunclegend->Draw("SAME");
        }
        
    }

    SystematicUncertainty_PtCanvas->SaveAs(Form("%s/SystematicUncertainty_Pt_Iter_%i.pdf", PlotDir.Data(),4));

    TCanvas *SystematicUncertainty_ZCanvas = new TCanvas("SystematicUncertainty_ZCanvas", "SystematicUncertainty_ZCanvas", 2100, 700);
    SystematicUncertainty_ZCanvas->Divide(5,3);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            SystematicUncertainty_ZCanvas->cd(cent*5+pT);
            for (int i = 0; i < nSystematic; i++){
                SystematicUncertainty_Z[pT-1][i][cent]->GetYaxis()->SetTitle(SystematicUncYaxisName.Data());
                SystematicUncertainty_Z[pT-1][i][cent]->GetYaxis()->SetRangeUser(-100,100);
                SystematicUncertainty_Z[pT-1][i][cent]->Draw("LP HIST SAME");
            }
            sysunclegend->Draw("SAME");
        }
        
    }

    SystematicUncertainty_ZCanvas->SaveAs(Form("%s/SystematicUncertainty_Z_Iter_%i.pdf", PlotDir.Data(),4));

    TCanvas *SystematicUncertainty_ZPerJetCanvas = new TCanvas("SystematicUncertainty_ZPerJetCanvas", "SystematicUncertainty_ZPerJetCanvas", 2100, 700);
    SystematicUncertainty_ZPerJetCanvas->Divide(5,3);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            SystematicUncertainty_ZPerJetCanvas->cd(cent*5+pT);
            for (int i = 0; i < nSystematic; i++){
                SystematicUncertainty_ZPerJet[pT-1][i][cent]->GetYaxis()->SetTitle(SystematicUncYaxisName.Data());
                SystematicUncertainty_ZPerJet[pT-1][i][cent]->GetYaxis()->SetRangeUser(-100,100);
                SystematicUncertainty_ZPerJet[pT-1][i][cent]->Draw("LP HIST SAME");
            }
            sysunclegend->Draw("SAME");
        }
        
    }

    SystematicUncertainty_ZPerJetCanvas->SaveAs(Form("%s/SystematicUncertainty_ZPerJet_Iter_%i.pdf", PlotDir.Data(),4));

    TCanvas *SystematicUncertainty_dRCanvas = new TCanvas("SystematicUncertainty_dRCanvas", "SystematicUncertainty_dRCanvas", 2100, 700);
    SystematicUncertainty_dRCanvas->Divide(5,3);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            SystematicUncertainty_dRCanvas->cd(cent*5+pT);
            for (int i = 0; i < nSystematic; i++){
                SystematicUncertainty_dR[pT-1][i][cent]->GetYaxis()->SetTitle(SystematicUncYaxisName.Data());
                SystematicUncertainty_dR[pT-1][i][cent]->GetYaxis()->SetRangeUser(-100,100);
                SystematicUncertainty_dR[pT-1][i][cent]->GetXaxis()->SetRangeUser(0,0.2);
                SystematicUncertainty_dR[pT-1][i][cent]->Draw("LP HIST SAME");
            }
            sysunclegend->Draw("SAME");
        }
        
    }

    SystematicUncertainty_dRCanvas->SaveAs(Form("%s/SystematicUncertainty_dR_Iter_%i.pdf", PlotDir.Data(),4));

    // ============== //

    TCanvas *SystematicUncertainty_RCP_PtCanvas = new TCanvas("SystematicUncertainty_RCP_PtCanvas", "SystematicUncertainty_RCP_PtCanvas", 2100, 700);
    SystematicUncertainty_RCP_PtCanvas->Divide(5,2);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){
            SystematicUncertainty_RCP_PtCanvas->cd(cent*5+pT);
            for (int i = 0; i < nSystematic; i++){
                SystematicUncertainty_RCP_Pt[pT-1][i][cent]->GetYaxis()->SetTitle(SystematicUncYaxisName.Data());
                SystematicUncertainty_RCP_Pt[pT-1][i][cent]->GetYaxis()->SetRangeUser(-100,100);
                SystematicUncertainty_RCP_Pt[pT-1][i][cent]->GetXaxis()->SetRangeUser(5,20);
                SystematicUncertainty_RCP_Pt[pT-1][i][cent]->Draw("LP HIST SAME");
                // if (cent == 0 && pT == 1) sysunclegend->AddEntry(SystematicUncertainty_RCP_Pt[pT-1][i][cent], SystematicName[i].Data(), "lp");
            }
            sysunclegend->Draw("SAME");
        }
        
    }

    SystematicUncertainty_RCP_PtCanvas->SaveAs(Form("%s/SystematicUncertainty_RCP_Pt_Iter_%i.pdf", PlotDir.Data(),4));

    TCanvas *SystematicUncertainty_RCP_ZCanvas = new TCanvas("SystematicUncertainty_RCP_ZCanvas", "SystematicUncertainty_RCP_ZCanvas", 2100, 700);
    SystematicUncertainty_RCP_ZCanvas->Divide(5,2);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){
            SystematicUncertainty_RCP_ZCanvas->cd(cent*5+pT);
            for (int i = 0; i < nSystematic; i++){
                SystematicUncertainty_RCP_Z[pT-1][i][cent]->GetYaxis()->SetTitle(SystematicUncYaxisName.Data());
                SystematicUncertainty_RCP_Z[pT-1][i][cent]->GetYaxis()->SetRangeUser(-100,100);
                SystematicUncertainty_RCP_Z[pT-1][i][cent]->Draw("LP HIST SAME");
            }
            sysunclegend->Draw("SAME");
        }
        
    }

    SystematicUncertainty_RCP_ZCanvas->SaveAs(Form("%s/SystematicUncertainty_RCP_Z_Iter_%i.pdf", PlotDir.Data(),4));

    TCanvas *SystematicUncertainty_RCP_ZPerJetCanvas = new TCanvas("SystematicUncertainty_RCP_ZPerJetCanvas", "SystematicUncertainty_RCP_ZPerJetCanvas", 2100, 700);
    SystematicUncertainty_RCP_ZPerJetCanvas->Divide(5,2);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){
            SystematicUncertainty_RCP_ZPerJetCanvas->cd(cent*5+pT);
            for (int i = 0; i < nSystematic; i++){
                SystematicUncertainty_RCP_ZPerJet[pT-1][i][cent]->GetYaxis()->SetTitle(SystematicUncYaxisName.Data());
                SystematicUncertainty_RCP_ZPerJet[pT-1][i][cent]->GetYaxis()->SetRangeUser(-100,100);
                SystematicUncertainty_RCP_ZPerJet[pT-1][i][cent]->Draw("LP HIST SAME");
            }
            sysunclegend->Draw("SAME");
        }
        
    }

    SystematicUncertainty_RCP_ZPerJetCanvas->SaveAs(Form("%s/SystematicUncertainty_RCP_ZPerJet_Iter_%i.pdf", PlotDir.Data(),4));

    TCanvas *SystematicUncertainty_RCP_dRCanvas = new TCanvas("SystematicUncertainty_RCP_dRCanvas", "SystematicUncertainty_RCP_dRCanvas", 2100, 700);
    SystematicUncertainty_RCP_dRCanvas->Divide(5,2);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){
            SystematicUncertainty_RCP_dRCanvas->cd(cent*5+pT);
            for (int i = 0; i < nSystematic; i++){
                SystematicUncertainty_RCP_dR[pT-1][i][cent]->GetYaxis()->SetTitle(SystematicUncYaxisName.Data());
                SystematicUncertainty_RCP_dR[pT-1][i][cent]->GetYaxis()->SetRangeUser(-100,100);
                SystematicUncertainty_RCP_dR[pT-1][i][cent]->GetXaxis()->SetRangeUser(0,0.2);
                SystematicUncertainty_RCP_dR[pT-1][i][cent]->Draw("LP HIST SAME");
            }
            sysunclegend->Draw("SAME");
        }
        
    }

    SystematicUncertainty_RCP_dRCanvas->SaveAs(Form("%s/SystematicUncertainty_RCP_dR_Iter_%i.pdf", PlotDir.Data(),4));

    // This part is the minimum and maximum systematic uncertainty variation for each bin from different iterations

    TH1D *MinSystematicUncertainty_Pt[5][3];
    TH1D *MinSystematicUncertainty_Z[5][3];
    TH1D *MinSystematicUncertainty_ZPerJet[5][3];
    TH1D *MinSystematicUncertainty_dR[5][3];

    TH1D *MaxSystematicUncertainty_Pt[5][3];
    TH1D *MaxSystematicUncertainty_Z[5][3];
    TH1D *MaxSystematicUncertainty_ZPerJet[5][3];
    TH1D *MaxSystematicUncertainty_dR[5][3];

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            MinSystematicUncertainty_Pt[pT-1][cent] = (TH1D *)JetpT[pT-1][0][cent]->Clone();
            MaxSystematicUncertainty_Pt[pT-1][cent] = (TH1D *)JetpT[pT-1][0][cent]->Clone();
            SetName(MinSystematicUncertainty_Pt[pT-1][cent], Form("%s Min Systematic Unc.", JetpT[pT-1][0][cent]->GetName()));
            SetName(MaxSystematicUncertainty_Pt[pT-1][cent], Form("%s Max Systematic Unc.", JetpT[pT-1][0][cent]->GetName()));
            for (int binx = 1; binx <= JetpT[pT-1][0][cent]->GetNbinsX(); binx++){
                double min = 0;
                double max = 0;
                for (int i = 0; i < 2; i++){
                    // cout << pT << "\t" << cent << "\t" << binx << "\t" << SystematicUncertainty_Pt[pT-1][i][cent]->GetBinContent(binx) << "\t" << min << "\t" << max << endl;
                    if (SystematicUncertainty_Pt[pT-1][i][cent]->GetBinContent(binx) < min) min = SystematicUncertainty_Pt[pT-1][i][cent]->GetBinContent(binx);
                    if (SystematicUncertainty_Pt[pT-1][i][cent]->GetBinContent(binx) > max) max = SystematicUncertainty_Pt[pT-1][i][cent]->GetBinContent(binx);
                }

                MaxSystematicUncertainty_Pt[pT-1][cent]->SetBinContent(binx, max);
                MaxSystematicUncertainty_Pt[pT-1][cent]->SetBinError(binx, 0);

                MinSystematicUncertainty_Pt[pT-1][cent]->SetBinContent(binx, min);
                MinSystematicUncertainty_Pt[pT-1][cent]->SetBinError(binx, 0);
            }
        }
    }

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            MinSystematicUncertainty_Z[pT-1][cent] = (TH1D *)JetZ[pT-1][0][cent]->Clone();
            MaxSystematicUncertainty_Z[pT-1][cent] = (TH1D *)JetZ[pT-1][0][cent]->Clone();
            SetName(MinSystematicUncertainty_Z[pT-1][cent], Form("%s Min Systematic Unc.", JetZ[pT-1][0][cent]->GetName()));
            SetName(MaxSystematicUncertainty_Z[pT-1][cent], Form("%s Max Systematic Unc.", JetZ[pT-1][0][cent]->GetName()));
            for (int binx = 1; binx <= JetZ[pT-1][0][cent]->GetNbinsX(); binx++){
                double min = 0;
                double max = 0;
                for (int i = 0; i < 2; i++){
                    if (SystematicUncertainty_Z[pT-1][i][cent]->GetBinContent(binx) < min) min = SystematicUncertainty_Z[pT-1][i][cent]->GetBinContent(binx);
                    if (SystematicUncertainty_Z[pT-1][i][cent]->GetBinContent(binx) > max) max = SystematicUncertainty_Z[pT-1][i][cent]->GetBinContent(binx);
                }

                MaxSystematicUncertainty_Z[pT-1][cent]->SetBinContent(binx, max);
                MaxSystematicUncertainty_Z[pT-1][cent]->SetBinError(binx, 0);

                MinSystematicUncertainty_Z[pT-1][cent]->SetBinContent(binx, min);
                MinSystematicUncertainty_Z[pT-1][cent]->SetBinError(binx, 0);

            }
        }
    }

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            MinSystematicUncertainty_ZPerJet[pT-1][cent] = (TH1D *)JetZPerJet[pT-1][0][cent]->Clone();
            MaxSystematicUncertainty_ZPerJet[pT-1][cent] = (TH1D *)JetZPerJet[pT-1][0][cent]->Clone();
            SetName(MinSystematicUncertainty_ZPerJet[pT-1][cent], Form("%s Min Systematic Unc.", JetZPerJet[pT-1][0][cent]->GetName()));
            SetName(MaxSystematicUncertainty_ZPerJet[pT-1][cent], Form("%s Max Systematic Unc.", JetZPerJet[pT-1][0][cent]->GetName()));
            for (int binx = 1; binx <= JetZPerJet[pT-1][0][cent]->GetNbinsX(); binx++){
                double min = 0;
                double max = 0;
                for (int i = 0; i < 2; i++){
                    if (SystematicUncertainty_ZPerJet[pT-1][i][cent]->GetBinContent(binx) < min) min = SystematicUncertainty_ZPerJet[pT-1][i][cent]->GetBinContent(binx);
                    if (SystematicUncertainty_ZPerJet[pT-1][i][cent]->GetBinContent(binx) > max) max = SystematicUncertainty_ZPerJet[pT-1][i][cent]->GetBinContent(binx);
                }

                MaxSystematicUncertainty_ZPerJet[pT-1][cent]->SetBinContent(binx, max);
                MaxSystematicUncertainty_ZPerJet[pT-1][cent]->SetBinError(binx, 0);

                MinSystematicUncertainty_ZPerJet[pT-1][cent]->SetBinContent(binx, min);
                MinSystematicUncertainty_ZPerJet[pT-1][cent]->SetBinError(binx, 0);

            }
        }
    }

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            MinSystematicUncertainty_dR[pT-1][cent] = (TH1D *)JetdR[pT-1][0][cent]->Clone();
            MaxSystematicUncertainty_dR[pT-1][cent] = (TH1D *)JetdR[pT-1][0][cent]->Clone();
            SetName(MinSystematicUncertainty_dR[pT-1][cent], Form("%s Min Systematic Unc.", JetdR[pT-1][0][cent]->GetName()));
            SetName(MaxSystematicUncertainty_dR[pT-1][cent], Form("%s Max Systematic Unc.", JetdR[pT-1][0][cent]->GetName()));
            for (int binx = 1; binx <= JetdR[pT-1][0][cent]->GetNbinsX(); binx++){
                double min = 0;
                double max = 0;
                for (int i = 0; i < 2; i++){
                    if (SystematicUncertainty_dR[pT-1][i][cent]->GetBinContent(binx) < min) min = SystematicUncertainty_dR[pT-1][i][cent]->GetBinContent(binx);
                    if (SystematicUncertainty_dR[pT-1][i][cent]->GetBinContent(binx) > max) max = SystematicUncertainty_dR[pT-1][i][cent]->GetBinContent(binx);
                }

                MaxSystematicUncertainty_dR[pT-1][cent]->SetBinContent(binx, max);
                MaxSystematicUncertainty_dR[pT-1][cent]->SetBinError(binx, 0);

                MinSystematicUncertainty_dR[pT-1][cent]->SetBinContent(binx, min);
                MinSystematicUncertainty_dR[pT-1][cent]->SetBinError(binx, 0);
            }
        }
    }

    auto minmaxiterlegend = new TLegend(0.1,0.1,0.8,0.3);

    TCanvas *MinMaxSystematicUncertainty_PtCanvas = new TCanvas("MinMaxSystematicUncertainty_PtCanvas", "MinMaxSystematicUncertainty_PtCanvas", 2100, 700);
    MinMaxSystematicUncertainty_PtCanvas->Divide(5,3);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            MinMaxSystematicUncertainty_PtCanvas->cd(cent*5+pT);
            // gPad->SetLogy();
            MaxSystematicUncertainty_Pt[pT-1][cent]->GetYaxis()->SetRangeUser(-100,100);
            MinSystematicUncertainty_Pt[pT-1][cent]->GetYaxis()->SetRangeUser(-100,100);
            MinSystematicUncertainty_Pt[pT-1][cent]->SetLineStyle(kDashed);
            MaxSystematicUncertainty_Pt[pT-1][cent]->GetYaxis()->SetTitle(SystematicUncYaxisName.Data());
            // MaxSystematicUncertainty_Pt[pT-1][cent]->GetXaxis()->SetRangeUser(5,20);
            // JetpT[pT-1][0][cent]->Draw("E2");
            MaxSystematicUncertainty_Pt[pT-1][cent]->Draw("E1");
            MinSystematicUncertainty_Pt[pT-1][cent]->Draw("E1 SAME");
            if (cent == 0 && pT == 1)minmaxiterlegend->AddEntry(MaxSystematicUncertainty_Pt[pT-1][cent], "Upper Limit Iter Systematic Unc.", "l");
            if (cent == 0 && pT == 1)minmaxiterlegend->AddEntry(MinSystematicUncertainty_Pt[pT-1][cent], "Lower Limit Iter Systematic Unc.", "l");

            minmaxiterlegend->Draw("SAME");
        }   
    }

    MinMaxSystematicUncertainty_PtCanvas->SaveAs(Form("%s/MinMaxSystematicUncertainty_Pt_IterOnly_Iter_%i.pdf", PlotDir.Data(),4));

    TCanvas *MinMaxSystematicUncertainty_ZCanvas = new TCanvas("MinMaxSystematicUncertainty_ZCanvas", "MinMaxSystematicUncertainty_ZCanvas", 2100, 700);
    MinMaxSystematicUncertainty_ZCanvas->Divide(5,3);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            MinMaxSystematicUncertainty_ZCanvas->cd(cent*5+pT);
            // gPad->SetLogy();
            MaxSystematicUncertainty_Z[pT-1][cent]->GetYaxis()->SetRangeUser(-100,100);
            MinSystematicUncertainty_Z[pT-1][cent]->GetYaxis()->SetRangeUser(-100,100);
            MinSystematicUncertainty_Z[pT-1][cent]->SetLineStyle(kDashed);
            MaxSystematicUncertainty_Z[pT-1][cent]->GetYaxis()->SetTitle(SystematicUncYaxisName.Data());
            // MaxSystematicUncertainty_Z[pT-1][cent]->GetXaxis()->SetRangeUser(5,20);
            // JetZ[pT-1][0][cent]->Draw("E2");
            MaxSystematicUncertainty_Z[pT-1][cent]->Draw("E1");
            MinSystematicUncertainty_Z[pT-1][cent]->Draw("E1 SAME");

            minmaxiterlegend->Draw("SAME");
        }   
    }

    MinMaxSystematicUncertainty_ZCanvas->SaveAs(Form("%s/MinMaxSystematicUncertainty_Z_IterOnly_Iter_%i.pdf", PlotDir.Data(),4));

    TCanvas *MinMaxSystematicUncertainty_ZPerJetCanvas = new TCanvas("MinMaxSystematicUncertainty_ZPerJetCanvas", "MinMaxSystematicUncertainty_ZPerJetCanvas", 2100, 700);
    MinMaxSystematicUncertainty_ZPerJetCanvas->Divide(5,3);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            MinMaxSystematicUncertainty_ZPerJetCanvas->cd(cent*5+pT);
            // gPad->SetLogy();
            MaxSystematicUncertainty_ZPerJet[pT-1][cent]->GetYaxis()->SetRangeUser(-100,100);
            MinSystematicUncertainty_ZPerJet[pT-1][cent]->GetYaxis()->SetRangeUser(-100,100);
            MinSystematicUncertainty_ZPerJet[pT-1][cent]->SetLineStyle(kDashed);
            MaxSystematicUncertainty_ZPerJet[pT-1][cent]->GetYaxis()->SetTitle(SystematicUncYaxisName.Data());
            // MaxSystematicUncertainty_Z[pT-1][cent]->GetXaxis()->SetRangeUser(5,20);
            // JetZ[pT-1][0][cent]->Draw("E2");
            MaxSystematicUncertainty_ZPerJet[pT-1][cent]->Draw("E1");
            MinSystematicUncertainty_ZPerJet[pT-1][cent]->Draw("E1 SAME");

            minmaxiterlegend->Draw("SAME");
        }   
    }

    MinMaxSystematicUncertainty_ZPerJetCanvas->SaveAs(Form("%s/MinMaxSystematicUncertainty_ZPerJet_IterOnly_Iter_%i.pdf", PlotDir.Data(),4));

    TCanvas *MinMaxSystematicUncertainty_dRCanvas = new TCanvas("MinMaxSystematicUncertainty_dRCanvas", "MinMaxSystematicUncertainty_dRCanvas", 2100, 700);
    MinMaxSystematicUncertainty_dRCanvas->Divide(5,3);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            MinMaxSystematicUncertainty_dRCanvas->cd(cent*5+pT);
            // gPad->SetLogy();
            MaxSystematicUncertainty_dR[pT-1][cent]->GetYaxis()->SetRangeUser(-100,100);
            MinSystematicUncertainty_dR[pT-1][cent]->GetYaxis()->SetRangeUser(-100,100);
            MinSystematicUncertainty_dR[pT-1][cent]->SetLineStyle(kDashed);
            MaxSystematicUncertainty_dR[pT-1][cent]->GetYaxis()->SetTitle(SystematicUncYaxisName.Data());
            // MaxSystematicUncertainty_dR[pT-1][cent]->GetXaxis()->SetRangeUser(5,20);
            // JetdR[pT-1][0][cent]->Draw("E2");
            MaxSystematicUncertainty_dR[pT-1][cent]->Draw("E1");
            MinSystematicUncertainty_dR[pT-1][cent]->Draw("E1 SAME");

            minmaxiterlegend->Draw("SAME");
        }   
    }

    MinMaxSystematicUncertainty_dRCanvas->SaveAs(Form("%s/MinMaxSystematicUncertainty_dR_IterOnly_Iter_%i.pdf", PlotDir.Data(),4));

    // This part is the minimum and maximum systematic uncertainty variation for each bin from different priors

    TH1D *MinSystematicUncertainty_Prior_Pt[5][3];
    TH1D *MinSystematicUncertainty_Prior_Z[5][3];
    TH1D *MinSystematicUncertainty_Prior_ZPerJet[5][3];
    TH1D *MinSystematicUncertainty_Prior_dR[5][3];

    TH1D *MaxSystematicUncertainty_Prior_Pt[5][3];
    TH1D *MaxSystematicUncertainty_Prior_Z[5][3];
    TH1D *MaxSystematicUncertainty_Prior_ZPerJet[5][3];
    TH1D *MaxSystematicUncertainty_Prior_dR[5][3];

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            MinSystematicUncertainty_Prior_Pt[pT-1][cent] = (TH1D *)JetpT[pT-1][0][cent]->Clone();
            MaxSystematicUncertainty_Prior_Pt[pT-1][cent] = (TH1D *)JetpT[pT-1][0][cent]->Clone();
            SetName(MinSystematicUncertainty_Prior_Pt[pT-1][cent], Form("%s Min Prior Systematic Unc.", JetpT[pT-1][0][cent]->GetName()));
            SetName(MaxSystematicUncertainty_Prior_Pt[pT-1][cent], Form("%s Max Prior Systematic Unc.", JetpT[pT-1][0][cent]->GetName()));
            for (int binx = 1; binx <= JetpT[pT-1][0][cent]->GetNbinsX(); binx++){
                double min = 0;
                double max = 0;
                for (int i = 2; i < 4; i++){
                    if (SystematicUncertainty_Pt[pT-1][i][cent]->GetBinContent(binx) < min) min = SystematicUncertainty_Pt[pT-1][i][cent]->GetBinContent(binx);
                    if (SystematicUncertainty_Pt[pT-1][i][cent]->GetBinContent(binx) > max) max = SystematicUncertainty_Pt[pT-1][i][cent]->GetBinContent(binx);
                }
                // cout << min << "\t" << max << endl;
                MaxSystematicUncertainty_Prior_Pt[pT-1][cent]->SetBinContent(binx, max);
                MaxSystematicUncertainty_Prior_Pt[pT-1][cent]->SetBinError(binx, 0);

                MinSystematicUncertainty_Prior_Pt[pT-1][cent]->SetBinContent(binx, min);
                MinSystematicUncertainty_Prior_Pt[pT-1][cent]->SetBinError(binx, 0);
            }
        }
    }

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            MinSystematicUncertainty_Prior_Z[pT-1][cent] = (TH1D *)JetZ[pT-1][0][cent]->Clone();
            MaxSystematicUncertainty_Prior_Z[pT-1][cent] = (TH1D *)JetZ[pT-1][0][cent]->Clone();
            SetName(MinSystematicUncertainty_Prior_Z[pT-1][cent], Form("%s Min Prior Systematic Unc.", JetZ[pT-1][0][cent]->GetName()));
            SetName(MaxSystematicUncertainty_Prior_Z[pT-1][cent], Form("%s Max Prior Systematic Unc.", JetZ[pT-1][0][cent]->GetName()));
            for (int binx = 1; binx <= JetZ[pT-1][0][cent]->GetNbinsX(); binx++){
                double min = 0;
                double max = 0;
                for (int i = 2; i < 4; i++){
                    if (SystematicUncertainty_Z[pT-1][i][cent]->GetBinContent(binx) < min) min = SystematicUncertainty_Z[pT-1][i][cent]->GetBinContent(binx);
                    if (SystematicUncertainty_Z[pT-1][i][cent]->GetBinContent(binx) > max) max = SystematicUncertainty_Z[pT-1][i][cent]->GetBinContent(binx);
                }

                MaxSystematicUncertainty_Prior_Z[pT-1][cent]->SetBinContent(binx, max);
                MaxSystematicUncertainty_Prior_Z[pT-1][cent]->SetBinError(binx, 0);

                MinSystematicUncertainty_Prior_Z[pT-1][cent]->SetBinContent(binx, min);
                MinSystematicUncertainty_Prior_Z[pT-1][cent]->SetBinError(binx, 0);
            }
        }
    }

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            MinSystematicUncertainty_Prior_ZPerJet[pT-1][cent] = (TH1D *)JetZPerJet[pT-1][0][cent]->Clone();
            MaxSystematicUncertainty_Prior_ZPerJet[pT-1][cent] = (TH1D *)JetZPerJet[pT-1][0][cent]->Clone();
            SetName(MinSystematicUncertainty_Prior_ZPerJet[pT-1][cent], Form("%s Min Prior Systematic Unc.", JetZPerJet[pT-1][0][cent]->GetName()));
            SetName(MaxSystematicUncertainty_Prior_ZPerJet[pT-1][cent], Form("%s Max Prior Systematic Unc.", JetZPerJet[pT-1][0][cent]->GetName()));
            for (int binx = 1; binx <= JetZPerJet[pT-1][0][cent]->GetNbinsX(); binx++){
                double min = 0;
                double max = 0;
                for (int i = 2; i < 4; i++){
                    if (SystematicUncertainty_ZPerJet[pT-1][i][cent]->GetBinContent(binx) < min) min = SystematicUncertainty_ZPerJet[pT-1][i][cent]->GetBinContent(binx);
                    if (SystematicUncertainty_ZPerJet[pT-1][i][cent]->GetBinContent(binx) > max) max = SystematicUncertainty_ZPerJet[pT-1][i][cent]->GetBinContent(binx);
                }

                MaxSystematicUncertainty_Prior_ZPerJet[pT-1][cent]->SetBinContent(binx, max);
                MaxSystematicUncertainty_Prior_ZPerJet[pT-1][cent]->SetBinError(binx, 0);

                MinSystematicUncertainty_Prior_ZPerJet[pT-1][cent]->SetBinContent(binx, min);
                MinSystematicUncertainty_Prior_ZPerJet[pT-1][cent]->SetBinError(binx, 0);
            }
        }
    }

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            MinSystematicUncertainty_Prior_dR[pT-1][cent] = (TH1D *)JetdR[pT-1][0][cent]->Clone();
            MaxSystematicUncertainty_Prior_dR[pT-1][cent] = (TH1D *)JetdR[pT-1][0][cent]->Clone();
            SetName(MinSystematicUncertainty_Prior_dR[pT-1][cent], Form("%s Min Prior Systematic Unc.", JetdR[pT-1][0][cent]->GetName()));
            SetName(MaxSystematicUncertainty_Prior_dR[pT-1][cent], Form("%s Max Prior Systematic Unc.", JetdR[pT-1][0][cent]->GetName()));
            for (int binx = 1; binx <= JetdR[pT-1][0][cent]->GetNbinsX(); binx++){
                double min = 0;
                double max = 0;
                for (int i = 2; i < 4; i++){
                    if (SystematicUncertainty_dR[pT-1][i][cent]->GetBinContent(binx) < min) min = SystematicUncertainty_dR[pT-1][i][cent]->GetBinContent(binx);
                    if (SystematicUncertainty_dR[pT-1][i][cent]->GetBinContent(binx) > max) max = SystematicUncertainty_dR[pT-1][i][cent]->GetBinContent(binx);
                }

                MaxSystematicUncertainty_Prior_dR[pT-1][cent]->SetBinContent(binx, max);
                MaxSystematicUncertainty_Prior_dR[pT-1][cent]->SetBinError(binx, 0);

                MinSystematicUncertainty_Prior_dR[pT-1][cent]->SetBinContent(binx, min);
                MinSystematicUncertainty_Prior_dR[pT-1][cent]->SetBinError(binx, 0);
            }
        }
    }

    auto minmaxpriorlegend = new TLegend(0.1,0.1,0.8,0.3);

    TCanvas *MinMaxSystematicUncertainty_Prior_PtCanvas = new TCanvas("MinMaxSystematicUncertainty_Prior_PtCanvas", "MinMaxSystematicUncertainty_Prior_PtCanvas", 2100, 700);
    MinMaxSystematicUncertainty_Prior_PtCanvas->Divide(5,3);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            MinMaxSystematicUncertainty_Prior_PtCanvas->cd(cent*5+pT);
            // gPad->SetLogy();
            MaxSystematicUncertainty_Prior_Pt[pT-1][cent]->GetYaxis()->SetRangeUser(-100,100);
            MinSystematicUncertainty_Prior_Pt[pT-1][cent]->GetYaxis()->SetRangeUser(-100,100);
            MinSystematicUncertainty_Prior_Pt[pT-1][cent]->SetLineStyle(kDashed);
            MaxSystematicUncertainty_Prior_Pt[pT-1][cent]->GetYaxis()->SetTitle(SystematicUncYaxisName.Data());
            // MaxSystematicUncertainty_Prior_Pt[pT-1][cent]->GetXaxis()->SetRangeUser(5,20);
            // JetpT[pT-1][0][cent]->Draw("E2");
            SetColor(MaxSystematicUncertainty_Prior_Pt[pT-1][cent], kRed, 0, 0, kRed);
            SetColor(MinSystematicUncertainty_Prior_Pt[pT-1][cent], kRed, 0, 0, kRed);
            MaxSystematicUncertainty_Prior_Pt[pT-1][cent]->Draw("E1");
            MinSystematicUncertainty_Prior_Pt[pT-1][cent]->Draw("E1 SAME");
            if (cent == 0 && pT == 1) minmaxpriorlegend->AddEntry(MaxSystematicUncertainty_Prior_Pt[pT-1][cent], "Upper Limit Prior Systematic Unc.", "l");
            if (cent == 0 && pT == 1) minmaxpriorlegend->AddEntry(MinSystematicUncertainty_Prior_Pt[pT-1][cent], "Lower Limit Prior Systematic Unc.", "l");
            minmaxpriorlegend->Draw("SAME");
        }   
    }

    MinMaxSystematicUncertainty_Prior_PtCanvas->SaveAs(Form("%s/MinMaxSystematicUncertainty_Prior_Pt_IterOnly_Iter_%i.pdf", PlotDir.Data(),4));

    TCanvas *MinMaxSystematicUncertainty_Prior_ZCanvas = new TCanvas("MinMaxSystematicUncertainty_Prior_ZCanvas", "MinMaxSystematicUncertainty_Prior_ZCanvas", 2100, 700);
    MinMaxSystematicUncertainty_Prior_ZCanvas->Divide(5,3);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            MinMaxSystematicUncertainty_Prior_ZCanvas->cd(cent*5+pT);
            // gPad->SetLogy();
            MaxSystematicUncertainty_Prior_Z[pT-1][cent]->GetYaxis()->SetRangeUser(-100,100);
            MinSystematicUncertainty_Prior_Z[pT-1][cent]->GetYaxis()->SetRangeUser(-100,100);
            MinSystematicUncertainty_Prior_Z[pT-1][cent]->SetLineStyle(kDashed);
            MaxSystematicUncertainty_Prior_Z[pT-1][cent]->GetYaxis()->SetTitle(SystematicUncYaxisName.Data());
            // MaxSystematicUncertainty_Prior_Z[pT-1][cent]->GetXaxis()->SetRangeUser(5,20);
            // JetZ[pT-1][0][cent]->Draw("E2");
            SetColor(MaxSystematicUncertainty_Prior_Z[pT-1][cent], kRed, 0, 0, kRed);
            SetColor(MinSystematicUncertainty_Prior_Z[pT-1][cent], kRed, 0, 0, kRed);
            MaxSystematicUncertainty_Prior_Z[pT-1][cent]->Draw("E1");
            MinSystematicUncertainty_Prior_Z[pT-1][cent]->Draw("E1 SAME");
            minmaxpriorlegend->Draw("SAME");
        }   
    }

    MinMaxSystematicUncertainty_Prior_ZCanvas->SaveAs(Form("%s/MinMaxSystematicUncertainty_Prior_Z_IterOnly_Iter_%i.pdf", PlotDir.Data(),4));

    TCanvas *MinMaxSystematicUncertainty_Prior_ZPerJetCanvas = new TCanvas("MinMaxSystematicUncertainty_Prior_ZPerJetCanvas", "MinMaxSystematicUncertainty_Prior_ZPerJetCanvas", 2100, 700);
    MinMaxSystematicUncertainty_Prior_ZPerJetCanvas->Divide(5,3);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            MinMaxSystematicUncertainty_Prior_ZPerJetCanvas->cd(cent*5+pT);
            // gPad->SetLogy();
            MaxSystematicUncertainty_Prior_ZPerJet[pT-1][cent]->GetYaxis()->SetRangeUser(-100,100);
            MinSystematicUncertainty_Prior_ZPerJet[pT-1][cent]->GetYaxis()->SetRangeUser(-100,100);
            MinSystematicUncertainty_Prior_ZPerJet[pT-1][cent]->SetLineStyle(kDashed);
            MaxSystematicUncertainty_Prior_ZPerJet[pT-1][cent]->GetYaxis()->SetTitle(SystematicUncYaxisName.Data());
            // MaxSystematicUncertainty_Prior_Z[pT-1][cent]->GetXaxis()->SetRangeUser(5,20);
            // JetZ[pT-1][0][cent]->Draw("E2");
            SetColor(MaxSystematicUncertainty_Prior_ZPerJet[pT-1][cent], kRed, 0, 0, kRed);
            SetColor(MinSystematicUncertainty_Prior_ZPerJet[pT-1][cent], kRed, 0, 0, kRed);
            MaxSystematicUncertainty_Prior_ZPerJet[pT-1][cent]->Draw("E1");
            MinSystematicUncertainty_Prior_ZPerJet[pT-1][cent]->Draw("E1 SAME");
            minmaxpriorlegend->Draw("SAME");
        }   
    }

    MinMaxSystematicUncertainty_Prior_ZPerJetCanvas->SaveAs(Form("%s/MinMaxSystematicUncertainty_Prior_ZPerJet_IterOnly_Iter_%i.pdf", PlotDir.Data(),4));

    TCanvas *MinMaxSystematicUncertainty_Prior_dRCanvas = new TCanvas("MinMaxSystematicUncertainty_Prior_dRCanvas", "MinMaxSystematicUncertainty_Prior_dRCanvas", 2100, 700);
    MinMaxSystematicUncertainty_Prior_dRCanvas->Divide(5,3);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            MinMaxSystematicUncertainty_Prior_dRCanvas->cd(cent*5+pT);
            // gPad->SetLogy();
            MaxSystematicUncertainty_Prior_dR[pT-1][cent]->GetYaxis()->SetRangeUser(-100,100);
            MinSystematicUncertainty_Prior_dR[pT-1][cent]->GetYaxis()->SetRangeUser(-100,100);
            MinSystematicUncertainty_Prior_dR[pT-1][cent]->SetLineStyle(kDashed);
            MaxSystematicUncertainty_Prior_dR[pT-1][cent]->GetYaxis()->SetTitle(SystematicUncYaxisName.Data());
            // MaxSystematicUncertainty_Prior_dR[pT-1][cent]->GetXaxis()->SetRangeUser(5,20);
            // JetdR[pT-1][0][cent]->Draw("E2");
            SetColor(MaxSystematicUncertainty_Prior_dR[pT-1][cent], kRed, 0, 0, kRed);
            SetColor(MinSystematicUncertainty_Prior_dR[pT-1][cent], kRed, 0, 0, kRed);
            MaxSystematicUncertainty_Prior_dR[pT-1][cent]->Draw("E1");
            MinSystematicUncertainty_Prior_dR[pT-1][cent]->Draw("E1 SAME");
            minmaxpriorlegend->Draw("SAME");
        }   
    }

    MinMaxSystematicUncertainty_Prior_dRCanvas->SaveAs(Form("%s/MinMaxSystematicUncertainty_Prior_dR_IterOnly_Iter_%i.pdf", PlotDir.Data(),4));

    // This part is the minimum and maximum systematic uncertainty variation for each bin from different iterations for RCPs

    TH1D *MinSystematicUncertainty_RCP_Pt[5][2];
    TH1D *MinSystematicUncertainty_RCP_Z[5][2];
    TH1D *MinSystematicUncertainty_RCP_ZPerJet[5][2];
    TH1D *MinSystematicUncertainty_RCP_dR[5][2];

    TH1D *MaxSystematicUncertainty_RCP_Pt[5][2];
    TH1D *MaxSystematicUncertainty_RCP_Z[5][2];
    TH1D *MaxSystematicUncertainty_RCP_ZPerJet[5][2];
    TH1D *MaxSystematicUncertainty_RCP_dR[5][2];

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){
            MinSystematicUncertainty_RCP_Pt[pT-1][cent] = (TH1D *)RCP_Pt[pT-1][0][cent]->Clone();
            MaxSystematicUncertainty_RCP_Pt[pT-1][cent] = (TH1D *)RCP_Pt[pT-1][0][cent]->Clone();
            SetName(MinSystematicUncertainty_RCP_Pt[pT-1][cent], Form("%s Min Systematic Unc.", RCP_Pt[pT-1][0][cent]->GetName()));
            SetName(MaxSystematicUncertainty_RCP_Pt[pT-1][cent], Form("%s Max Systematic Unc.", RCP_Pt[pT-1][0][cent]->GetName()));
            for (int binx = 1; binx <= RCP_Pt[pT-1][0][cent]->GetNbinsX(); binx++){
                double min = 0;
                double max = 0;
                for (int i = 0; i < 2; i++){
                    // cout << pT << "\t" << cent << "\t" << binx << "\t" << SystematicUncertainty_Pt[pT-1][i][cent]->GetBinContent(binx) << "\t" << min << "\t" << max << endl;
                    if (SystematicUncertainty_RCP_Pt[pT-1][i][cent]->GetBinContent(binx) < min) min = SystematicUncertainty_RCP_Pt[pT-1][i][cent]->GetBinContent(binx);
                    if (SystematicUncertainty_RCP_Pt[pT-1][i][cent]->GetBinContent(binx) > max) max = SystematicUncertainty_RCP_Pt[pT-1][i][cent]->GetBinContent(binx);
                }

                MaxSystematicUncertainty_RCP_Pt[pT-1][cent]->SetBinContent(binx, max);
                MaxSystematicUncertainty_RCP_Pt[pT-1][cent]->SetBinError(binx, 0);

                MinSystematicUncertainty_RCP_Pt[pT-1][cent]->SetBinContent(binx, min);
                MinSystematicUncertainty_RCP_Pt[pT-1][cent]->SetBinError(binx, 0);
            }
        }
    }

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){
            MinSystematicUncertainty_RCP_Z[pT-1][cent] = (TH1D *)RCP_Z[pT-1][0][cent]->Clone();
            MaxSystematicUncertainty_RCP_Z[pT-1][cent] = (TH1D *)RCP_Z[pT-1][0][cent]->Clone();
            SetName(MinSystematicUncertainty_RCP_Z[pT-1][cent], Form("%s Min Systematic Unc.", RCP_Z[pT-1][0][cent]->GetName()));
            SetName(MaxSystematicUncertainty_RCP_Z[pT-1][cent], Form("%s Max Systematic Unc.", RCP_Z[pT-1][0][cent]->GetName()));
            for (int binx = 1; binx <= RCP_Z[pT-1][0][cent]->GetNbinsX(); binx++){
                double min = 0;
                double max = 0;
                for (int i = 0; i < 2; i++){
                    if (SystematicUncertainty_RCP_Z[pT-1][i][cent]->GetBinContent(binx) < min) min = SystematicUncertainty_RCP_Z[pT-1][i][cent]->GetBinContent(binx);
                    if (SystematicUncertainty_RCP_Z[pT-1][i][cent]->GetBinContent(binx) > max) max = SystematicUncertainty_RCP_Z[pT-1][i][cent]->GetBinContent(binx);
                }

                MaxSystematicUncertainty_RCP_Z[pT-1][cent]->SetBinContent(binx, max);
                MaxSystematicUncertainty_RCP_Z[pT-1][cent]->SetBinError(binx, 0);

                MinSystematicUncertainty_RCP_Z[pT-1][cent]->SetBinContent(binx, min);
                MinSystematicUncertainty_RCP_Z[pT-1][cent]->SetBinError(binx, 0);

            }
        }
    }

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){
            MinSystematicUncertainty_RCP_ZPerJet[pT-1][cent] = (TH1D *)RCP_ZPerJet[pT-1][0][cent]->Clone();
            MaxSystematicUncertainty_RCP_ZPerJet[pT-1][cent] = (TH1D *)RCP_ZPerJet[pT-1][0][cent]->Clone();
            SetName(MinSystematicUncertainty_RCP_ZPerJet[pT-1][cent], Form("%s Min Systematic Unc.", RCP_ZPerJet[pT-1][0][cent]->GetName()));
            SetName(MaxSystematicUncertainty_RCP_ZPerJet[pT-1][cent], Form("%s Max Systematic Unc.", RCP_ZPerJet[pT-1][0][cent]->GetName()));
            for (int binx = 1; binx <= RCP_ZPerJet[pT-1][0][cent]->GetNbinsX(); binx++){
                double min = 0;
                double max = 0;
                for (int i = 0; i < 2; i++){
                    if (SystematicUncertainty_RCP_ZPerJet[pT-1][i][cent]->GetBinContent(binx) < min) min = SystematicUncertainty_RCP_ZPerJet[pT-1][i][cent]->GetBinContent(binx);
                    if (SystematicUncertainty_RCP_ZPerJet[pT-1][i][cent]->GetBinContent(binx) > max) max = SystematicUncertainty_RCP_ZPerJet[pT-1][i][cent]->GetBinContent(binx);
                }

                MaxSystematicUncertainty_RCP_ZPerJet[pT-1][cent]->SetBinContent(binx, max);
                MaxSystematicUncertainty_RCP_ZPerJet[pT-1][cent]->SetBinError(binx, 0);

                MinSystematicUncertainty_RCP_ZPerJet[pT-1][cent]->SetBinContent(binx, min);
                MinSystematicUncertainty_RCP_ZPerJet[pT-1][cent]->SetBinError(binx, 0);

            }
        }
    }

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){
            MinSystematicUncertainty_RCP_dR[pT-1][cent] = (TH1D *)RCP_dR[pT-1][0][cent]->Clone();
            MaxSystematicUncertainty_RCP_dR[pT-1][cent] = (TH1D *)RCP_dR[pT-1][0][cent]->Clone();
            SetName(MinSystematicUncertainty_RCP_dR[pT-1][cent], Form("%s Min Systematic Unc.", RCP_dR[pT-1][0][cent]->GetName()));
            SetName(MaxSystematicUncertainty_RCP_dR[pT-1][cent], Form("%s Max Systematic Unc.", RCP_dR[pT-1][0][cent]->GetName()));
            for (int binx = 1; binx <= RCP_dR[pT-1][0][cent]->GetNbinsX(); binx++){
                double min = 0;
                double max = 0;
                for (int i = 0; i < 2; i++){
                    if (SystematicUncertainty_RCP_dR[pT-1][i][cent]->GetBinContent(binx) < min) min = SystematicUncertainty_RCP_dR[pT-1][i][cent]->GetBinContent(binx);
                    if (SystematicUncertainty_RCP_dR[pT-1][i][cent]->GetBinContent(binx) > max) max = SystematicUncertainty_RCP_dR[pT-1][i][cent]->GetBinContent(binx);
                }

                MaxSystematicUncertainty_RCP_dR[pT-1][cent]->SetBinContent(binx, max);
                MaxSystematicUncertainty_RCP_dR[pT-1][cent]->SetBinError(binx, 0);

                MinSystematicUncertainty_RCP_dR[pT-1][cent]->SetBinContent(binx, min);
                MinSystematicUncertainty_RCP_dR[pT-1][cent]->SetBinError(binx, 0);
            }
        }
    }

    TCanvas *MinMaxSystematicUncertainty_RCP_PtCanvas = new TCanvas("MinMaxSystematicUncertainty_RCP_PtCanvas", "MinMaxSystematicUncertainty_RCP_PtCanvas", 2100, 700);
    MinMaxSystematicUncertainty_RCP_PtCanvas->Divide(5,2);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){
            MinMaxSystematicUncertainty_RCP_PtCanvas->cd(cent*5+pT);
            // gPad->SetLogy();
            MaxSystematicUncertainty_RCP_Pt[pT-1][cent]->GetYaxis()->SetRangeUser(-100,100);
            MinSystematicUncertainty_RCP_Pt[pT-1][cent]->GetYaxis()->SetRangeUser(-100,100);
            MinSystematicUncertainty_RCP_Pt[pT-1][cent]->SetLineStyle(kDashed);
            MaxSystematicUncertainty_RCP_Pt[pT-1][cent]->GetYaxis()->SetTitle(SystematicUncYaxisName.Data());
            // MaxSystematicUncertainty_Pt[pT-1][cent]->GetXaxis()->SetRangeUser(5,20);
            // JetpT[pT-1][0][cent]->Draw("E2");
            MaxSystematicUncertainty_RCP_Pt[pT-1][cent]->Draw("E1");
            MinSystematicUncertainty_RCP_Pt[pT-1][cent]->Draw("E1 SAME");
            // if (cent == 0 && pT == 1)minmaxiterlegend->AddEntry(MaxSystematicUncertainty_RCP_Pt[pT-1][cent], "Upper Limit Iter Systematic Unc.", "l");
            // if (cent == 0 && pT == 1)minmaxiterlegend->AddEntry(MinSystematicUncertainty_RCP_Pt[pT-1][cent], "Lower Limit Iter Systematic Unc.", "l");

            minmaxiterlegend->Draw("SAME");
        }   
    }

    MinMaxSystematicUncertainty_RCP_PtCanvas->SaveAs(Form("%s/MinMaxSystematicUncertainty_RCP_Pt_IterOnly_Iter_%i.pdf", PlotDir.Data(),4));

    TCanvas *MinMaxSystematicUncertainty_RCP_ZCanvas = new TCanvas("MinMaxSystematicUncertainty_RCP_ZCanvas", "MinMaxSystematicUncertainty_RCP_ZCanvas", 2100, 700);
    MinMaxSystematicUncertainty_RCP_ZCanvas->Divide(5,2);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){
            MinMaxSystematicUncertainty_RCP_ZCanvas->cd(cent*5+pT);
            // gPad->SetLogy();
            MaxSystematicUncertainty_RCP_Z[pT-1][cent]->GetYaxis()->SetRangeUser(-100,100);
            MinSystematicUncertainty_RCP_Z[pT-1][cent]->GetYaxis()->SetRangeUser(-100,100);
            MinSystematicUncertainty_RCP_Z[pT-1][cent]->SetLineStyle(kDashed);
            MaxSystematicUncertainty_RCP_Z[pT-1][cent]->GetYaxis()->SetTitle(SystematicUncYaxisName.Data());
            // MaxSystematicUncertainty_Z[pT-1][cent]->GetXaxis()->SetRangeUser(5,20);
            // JetZ[pT-1][0][cent]->Draw("E2");
            MaxSystematicUncertainty_RCP_Z[pT-1][cent]->Draw("E1");
            MinSystematicUncertainty_RCP_Z[pT-1][cent]->Draw("E1 SAME");

            minmaxiterlegend->Draw("SAME");
        }   
    }

    MinMaxSystematicUncertainty_RCP_ZCanvas->SaveAs(Form("%s/MinMaxSystematicUncertainty_RCP_Z_IterOnly_Iter_%i.pdf", PlotDir.Data(),4));

    TCanvas *MinMaxSystematicUncertainty_RCP_ZPerJetCanvas = new TCanvas("MinMaxSystematicUncertainty_RCP_ZPerJetCanvas", "MinMaxSystematicUncertainty_RCP_ZPerJetCanvas", 2100, 700);
    MinMaxSystematicUncertainty_RCP_ZPerJetCanvas->Divide(5,3);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            MinMaxSystematicUncertainty_RCP_ZPerJetCanvas->cd(cent*5+pT);
            // gPad->SetLogy();
            MaxSystematicUncertainty_RCP_ZPerJet[pT-1][cent]->GetYaxis()->SetRangeUser(-100,100);
            MinSystematicUncertainty_RCP_ZPerJet[pT-1][cent]->GetYaxis()->SetRangeUser(-100,100);
            MinSystematicUncertainty_RCP_ZPerJet[pT-1][cent]->SetLineStyle(kDashed);
            MaxSystematicUncertainty_RCP_ZPerJet[pT-1][cent]->GetYaxis()->SetTitle(SystematicUncYaxisName.Data());
            // MaxSystematicUncertainty_Z[pT-1][cent]->GetXaxis()->SetRangeUser(5,20);
            // JetZ[pT-1][0][cent]->Draw("E2");
            MaxSystematicUncertainty_RCP_ZPerJet[pT-1][cent]->Draw("E1");
            MinSystematicUncertainty_RCP_ZPerJet[pT-1][cent]->Draw("E1 SAME");

            minmaxiterlegend->Draw("SAME");
        }   
    }

    MinMaxSystematicUncertainty_RCP_ZPerJetCanvas->SaveAs(Form("%s/MinMaxSystematicUncertainty_RCP_ZPerJet_IterOnly_Iter_%i.pdf", PlotDir.Data(),4));

    TCanvas *MinMaxSystematicUncertainty_RCP_dRCanvas = new TCanvas("MinMaxSystematicUncertainty_RCP_dRCanvas", "MinMaxSystematicUncertainty_RCP_dRCanvas", 2100, 700);
    MinMaxSystematicUncertainty_RCP_dRCanvas->Divide(5,3);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            MinMaxSystematicUncertainty_RCP_dRCanvas->cd(cent*5+pT);
            // gPad->SetLogy();
            MaxSystematicUncertainty_RCP_dR[pT-1][cent]->GetYaxis()->SetRangeUser(-100,100);
            MinSystematicUncertainty_RCP_dR[pT-1][cent]->GetYaxis()->SetRangeUser(-100,100);
            MinSystematicUncertainty_RCP_dR[pT-1][cent]->SetLineStyle(kDashed);
            MaxSystematicUncertainty_RCP_dR[pT-1][cent]->GetYaxis()->SetTitle(SystematicUncYaxisName.Data());
            // MaxSystematicUncertainty_dR[pT-1][cent]->GetXaxis()->SetRangeUser(5,20);
            // JetdR[pT-1][0][cent]->Draw("E2");
            MaxSystematicUncertainty_RCP_dR[pT-1][cent]->Draw("E1");
            MinSystematicUncertainty_RCP_dR[pT-1][cent]->Draw("E1 SAME");

            minmaxiterlegend->Draw("SAME");
        }   
    }

    MinMaxSystematicUncertainty_RCP_dRCanvas->SaveAs(Form("%s/MinMaxSystematicUncertainty_RCP_dR_IterOnly_Iter_%i.pdf", PlotDir.Data(),4));

    // This part is the minimum and maximum systematic uncertainty variation for each bin from different priors

    TH1D *MinSystematicUncertainty_Prior_RCP_Pt[5][3];
    TH1D *MinSystematicUncertainty_Prior_RCP_Z[5][3];
    TH1D *MinSystematicUncertainty_Prior_RCP_ZPerJet[5][3];
    TH1D *MinSystematicUncertainty_Prior_RCP_dR[5][3];

    TH1D *MaxSystematicUncertainty_Prior_RCP_Pt[5][3];
    TH1D *MaxSystematicUncertainty_Prior_RCP_Z[5][3];
    TH1D *MaxSystematicUncertainty_Prior_RCP_ZPerJet[5][3];
    TH1D *MaxSystematicUncertainty_Prior_RCP_dR[5][3];

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){
            MinSystematicUncertainty_Prior_RCP_Pt[pT-1][cent] = (TH1D *)RCP_Pt[pT-1][0][cent]->Clone();
            MaxSystematicUncertainty_Prior_RCP_Pt[pT-1][cent] = (TH1D *)RCP_Pt[pT-1][0][cent]->Clone();
            SetName(MinSystematicUncertainty_Prior_RCP_Pt[pT-1][cent], Form("%s Min Prior Systematic Unc.", RCP_Pt[pT-1][0][cent]->GetName()));
            SetName(MaxSystematicUncertainty_Prior_RCP_Pt[pT-1][cent], Form("%s Max Prior Systematic Unc.", RCP_Pt[pT-1][0][cent]->GetName()));
            for (int binx = 1; binx <= RCP_Pt[pT-1][0][cent]->GetNbinsX(); binx++){
                double min = 0;
                double max = 0;
                for (int i = 2; i < 4; i++){
                    if (SystematicUncertainty_RCP_Pt[pT-1][i][cent]->GetBinContent(binx) < min) min = SystematicUncertainty_RCP_Pt[pT-1][i][cent]->GetBinContent(binx);
                    if (SystematicUncertainty_RCP_Pt[pT-1][i][cent]->GetBinContent(binx) > max) max = SystematicUncertainty_RCP_Pt[pT-1][i][cent]->GetBinContent(binx);
                }
                // cout << min << "\t" << max << endl;
                MaxSystematicUncertainty_Prior_RCP_Pt[pT-1][cent]->SetBinContent(binx, max);
                MaxSystematicUncertainty_Prior_RCP_Pt[pT-1][cent]->SetBinError(binx, 0);

                MinSystematicUncertainty_Prior_RCP_Pt[pT-1][cent]->SetBinContent(binx, min);
                MinSystematicUncertainty_Prior_RCP_Pt[pT-1][cent]->SetBinError(binx, 0);
            }
        }
    }

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){
            MinSystematicUncertainty_Prior_RCP_Z[pT-1][cent] = (TH1D *)RCP_Z[pT-1][0][cent]->Clone();
            MaxSystematicUncertainty_Prior_RCP_Z[pT-1][cent] = (TH1D *)RCP_Z[pT-1][0][cent]->Clone();
            SetName(MinSystematicUncertainty_Prior_RCP_Z[pT-1][cent], Form("%s Min Prior Systematic Unc.", RCP_Z[pT-1][0][cent]->GetName()));
            SetName(MaxSystematicUncertainty_Prior_RCP_Z[pT-1][cent], Form("%s Max Prior Systematic Unc.", RCP_Z[pT-1][0][cent]->GetName()));
            for (int binx = 1; binx <= RCP_Z[pT-1][0][cent]->GetNbinsX(); binx++){
                double min = 0;
                double max = 0;
                for (int i = 2; i < 4; i++){
                    if (SystematicUncertainty_RCP_Z[pT-1][i][cent]->GetBinContent(binx) < min) min = SystematicUncertainty_RCP_Z[pT-1][i][cent]->GetBinContent(binx);
                    if (SystematicUncertainty_RCP_Z[pT-1][i][cent]->GetBinContent(binx) > max) max = SystematicUncertainty_RCP_Z[pT-1][i][cent]->GetBinContent(binx);
                }

                MaxSystematicUncertainty_Prior_RCP_Z[pT-1][cent]->SetBinContent(binx, max);
                MaxSystematicUncertainty_Prior_RCP_Z[pT-1][cent]->SetBinError(binx, 0);

                MinSystematicUncertainty_Prior_RCP_Z[pT-1][cent]->SetBinContent(binx, min);
                MinSystematicUncertainty_Prior_RCP_Z[pT-1][cent]->SetBinError(binx, 0);
            }
        }
    }

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){
            MinSystematicUncertainty_Prior_RCP_ZPerJet[pT-1][cent] = (TH1D *)RCP_ZPerJet[pT-1][0][cent]->Clone();
            MaxSystematicUncertainty_Prior_RCP_ZPerJet[pT-1][cent] = (TH1D *)RCP_ZPerJet[pT-1][0][cent]->Clone();
            SetName(MinSystematicUncertainty_Prior_RCP_ZPerJet[pT-1][cent], Form("%s Min Prior Systematic Unc.", RCP_ZPerJet[pT-1][0][cent]->GetName()));
            SetName(MaxSystematicUncertainty_Prior_RCP_ZPerJet[pT-1][cent], Form("%s Max Prior Systematic Unc.", RCP_ZPerJet[pT-1][0][cent]->GetName()));
            for (int binx = 1; binx <= RCP_ZPerJet[pT-1][0][cent]->GetNbinsX(); binx++){
                double min = 0;
                double max = 0;
                for (int i = 2; i < 4; i++){
                    if (SystematicUncertainty_RCP_ZPerJet[pT-1][i][cent]->GetBinContent(binx) < min) min = SystematicUncertainty_RCP_ZPerJet[pT-1][i][cent]->GetBinContent(binx);
                    if (SystematicUncertainty_RCP_ZPerJet[pT-1][i][cent]->GetBinContent(binx) > max) max = SystematicUncertainty_RCP_ZPerJet[pT-1][i][cent]->GetBinContent(binx);
                }

                MaxSystematicUncertainty_Prior_RCP_ZPerJet[pT-1][cent]->SetBinContent(binx, max);
                MaxSystematicUncertainty_Prior_RCP_ZPerJet[pT-1][cent]->SetBinError(binx, 0);

                MinSystematicUncertainty_Prior_RCP_ZPerJet[pT-1][cent]->SetBinContent(binx, min);
                MinSystematicUncertainty_Prior_RCP_ZPerJet[pT-1][cent]->SetBinError(binx, 0);
            }
        }
    }

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){
            MinSystematicUncertainty_Prior_RCP_dR[pT-1][cent] = (TH1D *)RCP_dR[pT-1][0][cent]->Clone();
            MaxSystematicUncertainty_Prior_RCP_dR[pT-1][cent] = (TH1D *)RCP_dR[pT-1][0][cent]->Clone();
            SetName(MinSystematicUncertainty_Prior_RCP_dR[pT-1][cent], Form("%s Min Prior Systematic Unc.", RCP_dR[pT-1][0][cent]->GetName()));
            SetName(MaxSystematicUncertainty_Prior_RCP_dR[pT-1][cent], Form("%s Max Prior Systematic Unc.", RCP_dR[pT-1][0][cent]->GetName()));
            for (int binx = 1; binx <= RCP_dR[pT-1][0][cent]->GetNbinsX(); binx++){
                double min = 0;
                double max = 0;
                for (int i = 2; i < 4; i++){
                    if (SystematicUncertainty_RCP_dR[pT-1][i][cent]->GetBinContent(binx) < min) min = SystematicUncertainty_RCP_dR[pT-1][i][cent]->GetBinContent(binx);
                    if (SystematicUncertainty_RCP_dR[pT-1][i][cent]->GetBinContent(binx) > max) max = SystematicUncertainty_RCP_dR[pT-1][i][cent]->GetBinContent(binx);
                }

                MaxSystematicUncertainty_Prior_RCP_dR[pT-1][cent]->SetBinContent(binx, max);
                MaxSystematicUncertainty_Prior_RCP_dR[pT-1][cent]->SetBinError(binx, 0);

                MinSystematicUncertainty_Prior_RCP_dR[pT-1][cent]->SetBinContent(binx, min);
                MinSystematicUncertainty_Prior_RCP_dR[pT-1][cent]->SetBinError(binx, 0);
            }
        }
    }

    // auto minmaxpriorlegend = new TLegend(0.1,0.1,0.8,0.3);

    TCanvas *MinMaxSystematicUncertainty_Prior_RCP_PtCanvas = new TCanvas("MinMaxSystematicUncertainty_Prior_RCP_PtCanvas", "MinMaxSystematicUncertainty_Prior_RCP_PtCanvas", 2100, 700);
    MinMaxSystematicUncertainty_Prior_RCP_PtCanvas->Divide(5,2);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){
            MinMaxSystematicUncertainty_Prior_RCP_PtCanvas->cd(cent*5+pT);
            // gPad->SetLogy();
            MaxSystematicUncertainty_Prior_RCP_Pt[pT-1][cent]->GetYaxis()->SetRangeUser(-100,100);
            MinSystematicUncertainty_Prior_RCP_Pt[pT-1][cent]->GetYaxis()->SetRangeUser(-100,100);
            MinSystematicUncertainty_Prior_RCP_Pt[pT-1][cent]->SetLineStyle(kDashed);
            MaxSystematicUncertainty_Prior_RCP_Pt[pT-1][cent]->GetYaxis()->SetTitle(SystematicUncYaxisName.Data());
            // MaxSystematicUncertainty_Prior_Pt[pT-1][cent]->GetXaxis()->SetRangeUser(5,20);
            // JetpT[pT-1][0][cent]->Draw("E2");
            SetColor(MaxSystematicUncertainty_Prior_RCP_Pt[pT-1][cent], kRed, 0, 0, kRed);
            SetColor(MinSystematicUncertainty_Prior_RCP_Pt[pT-1][cent], kRed, 0, 0, kRed);
            MaxSystematicUncertainty_Prior_RCP_Pt[pT-1][cent]->Draw("E1");
            MinSystematicUncertainty_Prior_RCP_Pt[pT-1][cent]->Draw("E1 SAME");
            // if (cent == 0 && pT == 1) minmaxpriorlegend->AddEntry(MaxSystematicUncertainty_Prior_RCP_Pt[pT-1][cent], "Upper Limit Prior Systematic Unc.", "l");
            // if (cent == 0 && pT == 1) minmaxpriorlegend->AddEntry(MinSystematicUncertainty_Prior_RCP_Pt[pT-1][cent], "Lower Limit Prior Systematic Unc.", "l");
            minmaxpriorlegend->Draw("SAME");
        }   
    }

    MinMaxSystematicUncertainty_Prior_RCP_PtCanvas->SaveAs(Form("%s/MinMaxSystematicUncertainty_Prior_RCP_Pt_IterOnly_Iter_%i.pdf", PlotDir.Data(),4));

    TCanvas *MinMaxSystematicUncertainty_Prior_RCP_ZCanvas = new TCanvas("MinMaxSystematicUncertainty_Prior_RCP_ZCanvas", "MinMaxSystematicUncertainty_Prior_RCP_ZCanvas", 2100, 700);
    MinMaxSystematicUncertainty_Prior_RCP_ZCanvas->Divide(5,2);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){
            MinMaxSystematicUncertainty_Prior_RCP_ZCanvas->cd(cent*5+pT);
            // gPad->SetLogy();
            MaxSystematicUncertainty_Prior_RCP_Z[pT-1][cent]->GetYaxis()->SetRangeUser(-100,100);
            MinSystematicUncertainty_Prior_RCP_Z[pT-1][cent]->GetYaxis()->SetRangeUser(-100,100);
            MinSystematicUncertainty_Prior_RCP_Z[pT-1][cent]->SetLineStyle(kDashed);
            MaxSystematicUncertainty_Prior_RCP_Z[pT-1][cent]->GetYaxis()->SetTitle(SystematicUncYaxisName.Data());
            // MaxSystematicUncertainty_Prior_Z[pT-1][cent]->GetXaxis()->SetRangeUser(5,20);
            // JetZ[pT-1][0][cent]->Draw("E2");
            SetColor(MaxSystematicUncertainty_Prior_RCP_Z[pT-1][cent], kRed, 0, 0, kRed);
            SetColor(MinSystematicUncertainty_Prior_RCP_Z[pT-1][cent], kRed, 0, 0, kRed);
            MaxSystematicUncertainty_Prior_RCP_Z[pT-1][cent]->Draw("E1");
            MinSystematicUncertainty_Prior_RCP_Z[pT-1][cent]->Draw("E1 SAME");
            minmaxpriorlegend->Draw("SAME");
        }   
    }

    MinMaxSystematicUncertainty_Prior_RCP_ZCanvas->SaveAs(Form("%s/MinMaxSystematicUncertainty_Prior_RCP_Z_IterOnly_Iter_%i.pdf", PlotDir.Data(),4));

    TCanvas *MinMaxSystematicUncertainty_Prior_RCP_ZPerJetCanvas = new TCanvas("MinMaxSystematicUncertainty_Prior_RCP_ZPerJetCanvas", "MinMaxSystematicUncertainty_Prior_RCP_ZPerJetCanvas", 2100, 700);
    MinMaxSystematicUncertainty_Prior_RCP_ZPerJetCanvas->Divide(5,2);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){
            MinMaxSystematicUncertainty_Prior_RCP_ZPerJetCanvas->cd(cent*5+pT);
            // gPad->SetLogy();
            MaxSystematicUncertainty_Prior_RCP_ZPerJet[pT-1][cent]->GetYaxis()->SetRangeUser(-100,100);
            MinSystematicUncertainty_Prior_RCP_ZPerJet[pT-1][cent]->GetYaxis()->SetRangeUser(-100,100);
            MinSystematicUncertainty_Prior_RCP_ZPerJet[pT-1][cent]->SetLineStyle(kDashed);
            MaxSystematicUncertainty_Prior_RCP_ZPerJet[pT-1][cent]->GetYaxis()->SetTitle(SystematicUncYaxisName.Data());
            // MaxSystematicUncertainty_Prior_Z[pT-1][cent]->GetXaxis()->SetRangeUser(5,20);
            // JetZ[pT-1][0][cent]->Draw("E2");
            SetColor(MaxSystematicUncertainty_Prior_RCP_ZPerJet[pT-1][cent], kRed, 0, 0, kRed);
            SetColor(MinSystematicUncertainty_Prior_RCP_ZPerJet[pT-1][cent], kRed, 0, 0, kRed);
            MaxSystematicUncertainty_Prior_RCP_ZPerJet[pT-1][cent]->Draw("E1");
            MinSystematicUncertainty_Prior_RCP_ZPerJet[pT-1][cent]->Draw("E1 SAME");
            minmaxpriorlegend->Draw("SAME");
        }   
    }

    MinMaxSystematicUncertainty_Prior_RCP_ZPerJetCanvas->SaveAs(Form("%s/MinMaxSystematicUncertainty_Prior_RCP_ZPerJet_IterOnly_Iter_%i.pdf", PlotDir.Data(),4));

    TCanvas *MinMaxSystematicUncertainty_Prior_RCP_dRCanvas = new TCanvas("MinMaxSystematicUncertainty_Prior_RCP_dRCanvas", "MinMaxSystematicUncertainty_Prior_RCP_dRCanvas", 2100, 700);
    MinMaxSystematicUncertainty_Prior_RCP_dRCanvas->Divide(5,2);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){
            MinMaxSystematicUncertainty_Prior_RCP_dRCanvas->cd(cent*5+pT);
            // gPad->SetLogy();
            MaxSystematicUncertainty_Prior_RCP_dR[pT-1][cent]->GetYaxis()->SetRangeUser(-100,100);
            MinSystematicUncertainty_Prior_RCP_dR[pT-1][cent]->GetYaxis()->SetRangeUser(-100,100);
            MinSystematicUncertainty_Prior_RCP_dR[pT-1][cent]->SetLineStyle(kDashed);
            MaxSystematicUncertainty_Prior_RCP_dR[pT-1][cent]->GetYaxis()->SetTitle(SystematicUncYaxisName.Data());
            // MaxSystematicUncertainty_Prior_dR[pT-1][cent]->GetXaxis()->SetRangeUser(5,20);
            // JetdR[pT-1][0][cent]->Draw("E2");
            SetColor(MaxSystematicUncertainty_Prior_RCP_dR[pT-1][cent], kRed, 0, 0, kRed);
            SetColor(MinSystematicUncertainty_Prior_RCP_dR[pT-1][cent], kRed, 0, 0, kRed);
            MaxSystematicUncertainty_Prior_RCP_dR[pT-1][cent]->Draw("E1");
            MinSystematicUncertainty_Prior_RCP_dR[pT-1][cent]->Draw("E1 SAME");
            minmaxpriorlegend->Draw("SAME");
        }   
    }

    MinMaxSystematicUncertainty_Prior_RCP_dRCanvas->SaveAs(Form("%s/MinMaxSystematicUncertainty_Prior_RCP_dR_IterOnly_Iter_%i.pdf", PlotDir.Data(),4));

    // Finishing the minmax uncertainty for RCPs

    // Brining all the systematics together

    TH1D *Min_Iter_Pt[5][3];
    TH1D *Min_Iter_Z[5][3];
    TH1D *Min_Iter_ZPerJet[5][3];
    TH1D *Min_Iter_dR[5][3];

    TH1D *Max_Iter_Pt[5][3];
    TH1D *Max_Iter_Z[5][3];
    TH1D *Max_Iter_ZPerJet[5][3];
    TH1D *Max_Iter_dR[5][3];

    TH1D *Min_Prior_Pt[5][3];
    TH1D *Min_Prior_Z[5][3];
    TH1D *Min_Prior_ZPerJet[5][3];
    TH1D *Min_Prior_dR[5][3];

    TH1D *Max_Prior_Pt[5][3];
    TH1D *Max_Prior_Z[5][3];
    TH1D *Max_Prior_ZPerJet[5][3];
    TH1D *Max_Prior_dR[5][3];

    TH1D *Min_Track_Pt[5][3];
    TH1D *Min_Track_Z[5][3];
    TH1D *Min_Track_ZPerJet[5][3];
    TH1D *Min_Track_dR[5][3];

    TH1D *Max_Track_Pt[5][3];
    TH1D *Max_Track_Z[5][3];
    TH1D *Max_Track_ZPerJet[5][3];
    TH1D *Max_Track_dR[5][3];

    TH1D *Min_Total_Pt[5][3];
    TH1D *Min_Total_Z[5][3];
    TH1D *Min_Total_ZPerJet[5][3];
    TH1D *Min_Total_dR[5][3];

    TH1D *Max_Total_Pt[5][3];
    TH1D *Max_Total_Z[5][3];
    TH1D *Max_Total_ZPerJet[5][3];
    TH1D *Max_Total_dR[5][3];

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            Min_Iter_Pt[pT-1][cent] = (TH1D *)JetpT[pT-1][0][cent]->Clone();
            Max_Iter_Pt[pT-1][cent] = (TH1D *)JetpT[pT-1][0][cent]->Clone();
            SetName(Min_Iter_Pt[pT-1][cent], Form("%s Min With Iter Sys", JetpT[pT-1][0][cent]->GetName()));
            SetName(Max_Iter_Pt[pT-1][cent], Form("%s Max With Iter Sys", JetpT[pT-1][0][cent]->GetName()));

            Min_Prior_Pt[pT-1][cent] = (TH1D *)JetpT[pT-1][0][cent]->Clone();
            Max_Prior_Pt[pT-1][cent] = (TH1D *)JetpT[pT-1][0][cent]->Clone();
            SetName(Min_Prior_Pt[pT-1][cent], Form("%s Min With Prior Sys", JetpT[pT-1][0][cent]->GetName()));
            SetName(Max_Prior_Pt[pT-1][cent], Form("%s Max With Prior Sys", JetpT[pT-1][0][cent]->GetName()));

            Min_Track_Pt[pT-1][cent] = (TH1D *)JetpT[pT-1][0][cent]->Clone();
            Max_Track_Pt[pT-1][cent] = (TH1D *)JetpT[pT-1][0][cent]->Clone();
            SetName(Min_Track_Pt[pT-1][cent], Form("%s Min With Track Sys", JetpT[pT-1][0][cent]->GetName()));
            SetName(Max_Track_Pt[pT-1][cent], Form("%s Max With Track Sys", JetpT[pT-1][0][cent]->GetName()));

            Min_Total_Pt[pT-1][cent] = (TH1D *)JetpT[pT-1][0][cent]->Clone();
            Max_Total_Pt[pT-1][cent] = (TH1D *)JetpT[pT-1][0][cent]->Clone();
            SetName(Min_Total_Pt[pT-1][cent], Form("%s Min With Total Sys", JetpT[pT-1][0][cent]->GetName()));
            SetName(Max_Total_Pt[pT-1][cent], Form("%s Max With Total Sys", JetpT[pT-1][0][cent]->GetName()));

            for (int i = 1; i <= JetpT[pT-1][0][cent]->GetNbinsX(); i++){
                double minvaliter = JetpT[pT-1][0][cent]->GetBinContent(i) - MinSystematicUncertainty_Pt[pT-1][cent]->GetBinContent(i)*JetpT[pT-1][0][cent]->GetBinContent(i)/100;
                double maxvaliter = JetpT[pT-1][0][cent]->GetBinContent(i) + MaxSystematicUncertainty_Pt[pT-1][cent]->GetBinContent(i)*JetpT[pT-1][0][cent]->GetBinContent(i)/100;

                double minvalprior = JetpT[pT-1][0][cent]->GetBinContent(i) - MinSystematicUncertainty_Prior_Pt[pT-1][cent]->GetBinContent(i)*JetpT[pT-1][0][cent]->GetBinContent(i)/100;
                double maxvalprior = JetpT[pT-1][0][cent]->GetBinContent(i) + MaxSystematicUncertainty_Prior_Pt[pT-1][cent]->GetBinContent(i)*JetpT[pT-1][0][cent]->GetBinContent(i)/100;

                double minvaltotal = JetpT[pT-1][0][cent]->GetBinContent(i) - sqrt(pow(MinSystematicUncertainty_Pt[pT-1][cent]->GetBinContent(i)*JetpT[pT-1][0][cent]->GetBinContent(i)/100,2) + pow(MinSystematicUncertainty_Prior_Pt[pT-1][cent]->GetBinContent(i)*JetpT[pT-1][0][cent]->GetBinContent(i)/100,2));
                double maxvaltotal = JetpT[pT-1][0][cent]->GetBinContent(i) + sqrt(pow(MaxSystematicUncertainty_Pt[pT-1][cent]->GetBinContent(i)*JetpT[pT-1][0][cent]->GetBinContent(i)/100,2) + pow(MaxSystematicUncertainty_Prior_Pt[pT-1][cent]->GetBinContent(i)*JetpT[pT-1][0][cent]->GetBinContent(i)/100,2));

                Min_Iter_Pt[pT-1][cent]->SetBinContent(i, minvaliter);
                Min_Iter_Pt[pT-1][cent]->SetBinError(i, 0);

                Max_Iter_Pt[pT-1][cent]->SetBinContent(i, maxvaliter);
                Max_Iter_Pt[pT-1][cent]->SetBinError(i, 0);

                Min_Prior_Pt[pT-1][cent]->SetBinContent(i, minvalprior);
                Min_Prior_Pt[pT-1][cent]->SetBinError(i, 0);

                Max_Prior_Pt[pT-1][cent]->SetBinContent(i, maxvalprior);
                Max_Prior_Pt[pT-1][cent]->SetBinError(i, 0);

                Min_Total_Pt[pT-1][cent]->SetBinContent(i, minvaltotal);
                Min_Total_Pt[pT-1][cent]->SetBinError(i, 0);

                Max_Total_Pt[pT-1][cent]->SetBinContent(i, maxvaltotal);
                Max_Total_Pt[pT-1][cent]->SetBinError(i, 0);
            }
        }
    }

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            Min_Iter_Z[pT-1][cent] = (TH1D *)JetZ[pT-1][0][cent]->Clone();
            Max_Iter_Z[pT-1][cent] = (TH1D *)JetZ[pT-1][0][cent]->Clone();
            SetName(Min_Iter_Z[pT-1][cent], Form("%s Min With Iter Sys", JetZ[pT-1][0][cent]->GetName()));
            SetName(Max_Iter_Z[pT-1][cent], Form("%s Max With Iter Sys", JetZ[pT-1][0][cent]->GetName()));

            Min_Prior_Z[pT-1][cent] = (TH1D *)JetZ[pT-1][0][cent]->Clone();
            Max_Prior_Z[pT-1][cent] = (TH1D *)JetZ[pT-1][0][cent]->Clone();
            SetName(Min_Prior_Z[pT-1][cent], Form("%s Min With Prior Sys", JetZ[pT-1][0][cent]->GetName()));
            SetName(Max_Prior_Z[pT-1][cent], Form("%s Max With Prior Sys", JetZ[pT-1][0][cent]->GetName()));

            Min_Total_Z[pT-1][cent] = (TH1D *)JetZ[pT-1][0][cent]->Clone();
            Max_Total_Z[pT-1][cent] = (TH1D *)JetZ[pT-1][0][cent]->Clone();
            SetName(Min_Total_Z[pT-1][cent], Form("%s Min With Total Sys", JetZ[pT-1][0][cent]->GetName()));
            SetName(Max_Total_Z[pT-1][cent], Form("%s Max With Total Sys", JetZ[pT-1][0][cent]->GetName()));

            for (int i = 1; i <= JetZ[pT-1][0][cent]->GetNbinsX(); i++){
                double minvaliter = JetZ[pT-1][0][cent]->GetBinContent(i) - MinSystematicUncertainty_Z[pT-1][cent]->GetBinContent(i)*JetZ[pT-1][0][cent]->GetBinContent(i)/100;
                double maxvaliter = JetZ[pT-1][0][cent]->GetBinContent(i) + MaxSystematicUncertainty_Z[pT-1][cent]->GetBinContent(i)*JetZ[pT-1][0][cent]->GetBinContent(i)/100;

                double minvalprior = JetZ[pT-1][0][cent]->GetBinContent(i) - MinSystematicUncertainty_Prior_Z[pT-1][cent]->GetBinContent(i)*JetZ[pT-1][0][cent]->GetBinContent(i)/100;
                double maxvalprior = JetZ[pT-1][0][cent]->GetBinContent(i) + MaxSystematicUncertainty_Prior_Z[pT-1][cent]->GetBinContent(i)*JetZ[pT-1][0][cent]->GetBinContent(i)/100;

                double minvaltotal = JetZ[pT-1][0][cent]->GetBinContent(i) - sqrt(pow(MinSystematicUncertainty_Z[pT-1][cent]->GetBinContent(i)*JetZ[pT-1][0][cent]->GetBinContent(i)/100,2) + pow(MinSystematicUncertainty_Prior_Z[pT-1][cent]->GetBinContent(i)*JetZ[pT-1][0][cent]->GetBinContent(i)/100,2));
                double maxvaltotal = JetZ[pT-1][0][cent]->GetBinContent(i) + sqrt(pow(MaxSystematicUncertainty_Z[pT-1][cent]->GetBinContent(i)*JetZ[pT-1][0][cent]->GetBinContent(i)/100,2) + pow(MaxSystematicUncertainty_Prior_Z[pT-1][cent]->GetBinContent(i)*JetZ[pT-1][0][cent]->GetBinContent(i)/100,2));

                Min_Iter_Z[pT-1][cent]->SetBinContent(i, minvaliter);
                Min_Iter_Z[pT-1][cent]->SetBinError(i, 0);

                Max_Iter_Z[pT-1][cent]->SetBinContent(i, maxvaliter);
                Max_Iter_Z[pT-1][cent]->SetBinError(i, 0);

                Min_Prior_Z[pT-1][cent]->SetBinContent(i, minvalprior);
                Min_Prior_Z[pT-1][cent]->SetBinError(i, 0);

                Max_Prior_Z[pT-1][cent]->SetBinContent(i, maxvalprior);
                Max_Prior_Z[pT-1][cent]->SetBinError(i, 0);

                Min_Total_Z[pT-1][cent]->SetBinContent(i, minvaltotal);
                Min_Total_Z[pT-1][cent]->SetBinError(i, 0);

                Max_Total_Z[pT-1][cent]->SetBinContent(i, maxvaltotal);
                Max_Total_Z[pT-1][cent]->SetBinError(i, 0);
            }
        }
    }

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            Min_Iter_ZPerJet[pT-1][cent] = (TH1D *)JetZPerJet[pT-1][0][cent]->Clone();
            Max_Iter_ZPerJet[pT-1][cent] = (TH1D *)JetZPerJet[pT-1][0][cent]->Clone();
            SetName(Min_Iter_ZPerJet[pT-1][cent], Form("%s Min With Iter Sys", JetZPerJet[pT-1][0][cent]->GetName()));
            SetName(Max_Iter_ZPerJet[pT-1][cent], Form("%s Max With Iter Sys", JetZPerJet[pT-1][0][cent]->GetName()));

            Min_Prior_ZPerJet[pT-1][cent] = (TH1D *)JetZPerJet[pT-1][0][cent]->Clone();
            Max_Prior_ZPerJet[pT-1][cent] = (TH1D *)JetZPerJet[pT-1][0][cent]->Clone();
            SetName(Min_Prior_ZPerJet[pT-1][cent], Form("%s Min With Prior Sys", JetZPerJet[pT-1][0][cent]->GetName()));
            SetName(Max_Prior_ZPerJet[pT-1][cent], Form("%s Max With Prior Sys", JetZPerJet[pT-1][0][cent]->GetName()));

            Min_Total_ZPerJet[pT-1][cent] = (TH1D *)JetZPerJet[pT-1][0][cent]->Clone();
            Max_Total_ZPerJet[pT-1][cent] = (TH1D *)JetZPerJet[pT-1][0][cent]->Clone();
            SetName(Min_Total_ZPerJet[pT-1][cent], Form("%s Min With Total Sys", JetZPerJet[pT-1][0][cent]->GetName()));
            SetName(Max_Total_ZPerJet[pT-1][cent], Form("%s Max With Total Sys", JetZPerJet[pT-1][0][cent]->GetName()));

            for (int i = 1; i <= JetZPerJet[pT-1][0][cent]->GetNbinsX(); i++){
                double minvaliter = JetZPerJet[pT-1][0][cent]->GetBinContent(i) - MinSystematicUncertainty_ZPerJet[pT-1][cent]->GetBinContent(i)*JetZPerJet[pT-1][0][cent]->GetBinContent(i)/100;
                double maxvaliter = JetZPerJet[pT-1][0][cent]->GetBinContent(i) + MaxSystematicUncertainty_ZPerJet[pT-1][cent]->GetBinContent(i)*JetZPerJet[pT-1][0][cent]->GetBinContent(i)/100;

                double minvalprior = JetZPerJet[pT-1][0][cent]->GetBinContent(i) - MinSystematicUncertainty_Prior_ZPerJet[pT-1][cent]->GetBinContent(i)*JetZPerJet[pT-1][0][cent]->GetBinContent(i)/100;
                double maxvalprior = JetZPerJet[pT-1][0][cent]->GetBinContent(i) + MaxSystematicUncertainty_Prior_ZPerJet[pT-1][cent]->GetBinContent(i)*JetZPerJet[pT-1][0][cent]->GetBinContent(i)/100;

                double minvaltotal = JetZPerJet[pT-1][0][cent]->GetBinContent(i) - sqrt(pow(MinSystematicUncertainty_ZPerJet[pT-1][cent]->GetBinContent(i)*JetZPerJet[pT-1][0][cent]->GetBinContent(i)/100,2) + pow(MinSystematicUncertainty_Prior_ZPerJet[pT-1][cent]->GetBinContent(i)*JetZPerJet[pT-1][0][cent]->GetBinContent(i)/100,2));
                double maxvaltotal = JetZPerJet[pT-1][0][cent]->GetBinContent(i) + sqrt(pow(MaxSystematicUncertainty_ZPerJet[pT-1][cent]->GetBinContent(i)*JetZPerJet[pT-1][0][cent]->GetBinContent(i)/100,2) + pow(MaxSystematicUncertainty_Prior_ZPerJet[pT-1][cent]->GetBinContent(i)*JetZPerJet[pT-1][0][cent]->GetBinContent(i)/100,2));

                Min_Iter_ZPerJet[pT-1][cent]->SetBinContent(i, minvaliter);
                Min_Iter_ZPerJet[pT-1][cent]->SetBinError(i, 0);

                Max_Iter_ZPerJet[pT-1][cent]->SetBinContent(i, maxvaliter);
                Max_Iter_ZPerJet[pT-1][cent]->SetBinError(i, 0);

                Min_Prior_ZPerJet[pT-1][cent]->SetBinContent(i, minvalprior);
                Min_Prior_ZPerJet[pT-1][cent]->SetBinError(i, 0);

                Max_Prior_ZPerJet[pT-1][cent]->SetBinContent(i, maxvalprior);
                Max_Prior_ZPerJet[pT-1][cent]->SetBinError(i, 0);

                Min_Total_ZPerJet[pT-1][cent]->SetBinContent(i, minvaltotal);
                Min_Total_ZPerJet[pT-1][cent]->SetBinError(i, 0);

                Max_Total_ZPerJet[pT-1][cent]->SetBinContent(i, maxvaltotal);
                Max_Total_ZPerJet[pT-1][cent]->SetBinError(i, 0);
            }
        }
    }

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            Min_Iter_dR[pT-1][cent] = (TH1D *)JetdR[pT-1][0][cent]->Clone();
            Max_Iter_dR[pT-1][cent] = (TH1D *)JetdR[pT-1][0][cent]->Clone();
            SetName(Min_Iter_dR[pT-1][cent], Form("%s Min With Iter Sys", JetdR[pT-1][0][cent]->GetName()));
            SetName(Max_Iter_dR[pT-1][cent], Form("%s Max With Iter Sys", JetdR[pT-1][0][cent]->GetName()));

            Min_Prior_dR[pT-1][cent] = (TH1D *)JetdR[pT-1][0][cent]->Clone();
            Max_Prior_dR[pT-1][cent] = (TH1D *)JetdR[pT-1][0][cent]->Clone();
            SetName(Min_Prior_dR[pT-1][cent], Form("%s Min With Prior Sys", JetdR[pT-1][0][cent]->GetName()));
            SetName(Max_Prior_dR[pT-1][cent], Form("%s Max With Prior Sys", JetdR[pT-1][0][cent]->GetName()));

            Min_Total_dR[pT-1][cent] = (TH1D *)JetdR[pT-1][0][cent]->Clone();
            Max_Total_dR[pT-1][cent] = (TH1D *)JetdR[pT-1][0][cent]->Clone();
            SetName(Min_Total_dR[pT-1][cent], Form("%s Min With Total Sys", JetdR[pT-1][0][cent]->GetName()));
            SetName(Max_Total_dR[pT-1][cent], Form("%s Max With Total Sys", JetdR[pT-1][0][cent]->GetName()));

            for (int i = 1; i <= JetdR[pT-1][0][cent]->GetNbinsX(); i++){
                double minvaliter = JetdR[pT-1][0][cent]->GetBinContent(i) - MinSystematicUncertainty_dR[pT-1][cent]->GetBinContent(i)*JetdR[pT-1][0][cent]->GetBinContent(i)/100;
                double maxvaliter = JetdR[pT-1][0][cent]->GetBinContent(i) + MaxSystematicUncertainty_dR[pT-1][cent]->GetBinContent(i)*JetdR[pT-1][0][cent]->GetBinContent(i)/100;

                double minvalprior = JetdR[pT-1][0][cent]->GetBinContent(i) - MinSystematicUncertainty_Prior_dR[pT-1][cent]->GetBinContent(i)*JetdR[pT-1][0][cent]->GetBinContent(i)/100;
                double maxvalprior = JetdR[pT-1][0][cent]->GetBinContent(i) + MaxSystematicUncertainty_Prior_dR[pT-1][cent]->GetBinContent(i)*JetdR[pT-1][0][cent]->GetBinContent(i)/100;

                double minvaltotal = JetdR[pT-1][0][cent]->GetBinContent(i) - sqrt(pow(MinSystematicUncertainty_dR[pT-1][cent]->GetBinContent(i)*JetdR[pT-1][0][cent]->GetBinContent(i)/100,2) + pow(MinSystematicUncertainty_Prior_dR[pT-1][cent]->GetBinContent(i)*JetdR[pT-1][0][cent]->GetBinContent(i)/100,2));
                double maxvaltotal = JetdR[pT-1][0][cent]->GetBinContent(i) + sqrt(pow(MaxSystematicUncertainty_dR[pT-1][cent]->GetBinContent(i)*JetdR[pT-1][0][cent]->GetBinContent(i)/100,2) + pow(MaxSystematicUncertainty_Prior_dR[pT-1][cent]->GetBinContent(i)*JetdR[pT-1][0][cent]->GetBinContent(i)/100,2));

                Min_Iter_dR[pT-1][cent]->SetBinContent(i, minvaliter);
                Min_Iter_dR[pT-1][cent]->SetBinError(i, 0);

                Max_Iter_dR[pT-1][cent]->SetBinContent(i, maxvaliter);
                Max_Iter_dR[pT-1][cent]->SetBinError(i, 0);

                Min_Prior_dR[pT-1][cent]->SetBinContent(i, minvalprior);
                Min_Prior_dR[pT-1][cent]->SetBinError(i, 0);

                Max_Prior_dR[pT-1][cent]->SetBinContent(i, maxvalprior);
                Max_Prior_dR[pT-1][cent]->SetBinError(i, 0);

                Min_Total_dR[pT-1][cent]->SetBinContent(i, minvaltotal);
                Min_Total_dR[pT-1][cent]->SetBinError(i, 0);

                Max_Total_dR[pT-1][cent]->SetBinContent(i, maxvaltotal);
                Max_Total_dR[pT-1][cent]->SetBinError(i, 0);
            }
        }
    }



    TCanvas *MinMax_PtCanvas = new TCanvas("MinMax_Iter_PtCanvas", "MinMax_Iter_PtCanvas", 2100, 700);
    MinMax_PtCanvas->Divide(5,3);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            MinMax_PtCanvas->cd(cent*5+pT);
            gPad->SetLogy();
            Max_Total_Pt[pT-1][cent]->GetXaxis()->SetRangeUser(5,20);
            Max_Total_Pt[pT-1][cent]->GetYaxis()->SetRangeUser(1e-10,1e-3);

            SetColor(Max_Total_Pt[pT-1][cent], 0);
            Max_Total_Pt[pT-1][cent]->SetFillStyle(3344);
            Max_Total_Pt[pT-1][cent]->SetFillColorAlpha(kGreen-2, 0.6);
            Min_Total_Pt[pT-1][cent]->SetFillStyle(1001);
            Min_Total_Pt[pT-1][cent]->SetFillColor(10);
            Max_Total_Pt[pT-1][cent]->Draw("HIST ][");
            Min_Total_Pt[pT-1][cent]->Draw("SAME ][");

            AvgJetpT[pT-1][cent]->SetMarkerStyle(24);
            AvgJetpT[pT-1][cent]->Draw("E1 SAME");
            gPad->RedrawAxis();
        }
    }

    MinMax_PtCanvas->SaveAs(Form("%s/MinMax_Pt_Iter_%i.pdf", PlotDir.Data(),4));

    TFile *OldQMJetPtResults = new TFile("/Volumes/WorkDrive/work/2022/PreliminaryPlots/April2/JetSpectra.root", "READ");
    TH1D *JetPt_Central = (TH1D *)OldQMJetPtResults->Get("JetSpectra_0_10_Sys");
    TH1D *JetPt_MidCentral = (TH1D *)OldQMJetPtResults->Get("JetSpectra_10_40_Sys");
    TH1D *JetPt_Peripheral = (TH1D *)OldQMJetPtResults->Get("JetSpectra_40_80_Sys");

    auto QMJetPtlegend = new TLegend(0.6, 0.6, 0.9, 0.9);
    QMJetPtlegend->AddEntry(JetPt_Central, "QM 2022 Version", "lp");
    QMJetPtlegend->AddEntry(Max_Total_dR[4][0], "Current Version", "f");

    TCanvas *MinMax_PtCanvas_WithOldQMResults = (TCanvas *)MinMax_PtCanvas->DrawClone();

    MinMax_PtCanvas_WithOldQMResults->cd(5);
    SetColor(JetPt_Central, kRed, 20);
    JetPt_Central->Draw("EP SAME");
    QMJetPtlegend->Draw("SAME");
    MinMax_PtCanvas_WithOldQMResults->cd(10);
    JetPt_MidCentral->Draw("EP SAME");
    SetColor(JetPt_MidCentral, kRed, 20);
    QMJetPtlegend->Draw("SAME");
    MinMax_PtCanvas_WithOldQMResults->cd(15);
    SetColor(JetPt_Peripheral, kRed, 20);
    JetPt_Peripheral->Draw("EP SAME");
    QMJetPtlegend->Draw("SAME");

    MinMax_PtCanvas_WithOldQMResults->SaveAs(Form("%s/MinMax_Pt_Iter_%i_QM_Comparison.pdf", PlotDir.Data(),4));

    TCanvas *MinMax_ZCanvas = new TCanvas("MinMax_Iter_ZCanvas", "MinMax_Iter_ZCanvas", 2100, 700);
    MinMax_ZCanvas->Divide(5,3);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            MinMax_ZCanvas->cd(cent*5+pT);
            gPad->SetLogy();
            Max_Iter_Z[pT-1][cent]->GetXaxis()->SetRangeUser(5,20);
            Max_Iter_Z[pT-1][cent]->GetYaxis()->SetRangeUser(1e-10,1e-4);

            SetColor(Max_Total_Z[pT-1][cent], 0);
            Max_Total_Z[pT-1][cent]->SetFillStyle(3344);
            Max_Total_Z[pT-1][cent]->SetFillColorAlpha(kGreen-2, 0.6);
            Min_Total_Z[pT-1][cent]->SetFillStyle(1001);
            Min_Total_Z[pT-1][cent]->SetFillColor(10);
            Max_Total_Z[pT-1][cent]->Draw("HIST ][");
            Min_Total_Z[pT-1][cent]->Draw("SAME ][");

            AvgJetZ[pT-1][cent]->SetMarkerStyle(24);
            AvgJetZ[pT-1][cent]->Draw("E1 SAME");
            gPad->RedrawAxis();
        }
    }

    MinMax_ZCanvas->SaveAs(Form("%s/MinMax_Z_Iter_%i.pdf", PlotDir.Data(),4));

    TCanvas *MinMax_ZPerJetCanvas = new TCanvas("MinMax_Iter_ZPerJetCanvas", "MinMax_Iter_ZPerJetCanvas", 2100, 700);
    MinMax_ZPerJetCanvas->Divide(5,3);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            MinMax_ZPerJetCanvas->cd(cent*5+pT);
            gPad->SetLogy();
            Max_Iter_ZPerJet[pT-1][cent]->GetXaxis()->SetRangeUser(5,20);
            Max_Iter_ZPerJet[pT-1][cent]->GetYaxis()->SetRangeUser(1e-4,1e4);

            SetColor(Max_Total_ZPerJet[pT-1][cent], 0);
            Max_Total_ZPerJet[pT-1][cent]->SetFillStyle(3344);
            Max_Total_ZPerJet[pT-1][cent]->SetFillColorAlpha(kGreen-2, 0.6);
            Min_Total_ZPerJet[pT-1][cent]->SetFillStyle(1001);
            Min_Total_ZPerJet[pT-1][cent]->SetFillColor(10);
            Max_Total_ZPerJet[pT-1][cent]->Draw("HIST ][");
            Min_Total_ZPerJet[pT-1][cent]->Draw("SAME ][");

            AvgJetZPerJet[pT-1][cent]->SetMarkerStyle(24);
            AvgJetZPerJet[pT-1][cent]->Draw("E1 SAME");
            gPad->RedrawAxis();
        }
    }

    MinMax_ZPerJetCanvas->SaveAs(Form("%s/MinMax_ZPerJet_Iter_%i.pdf", PlotDir.Data(),4));

    TCanvas *MinMax_dRCanvas = new TCanvas("MinMax_dRCanvas", "MinMax_dRCanvas", 2100, 700);
    MinMax_dRCanvas->Divide(5,3);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            MinMax_dRCanvas->cd(cent*5+pT);
            gPad->SetLogy();
            Max_Iter_dR[pT-1][cent]->GetXaxis()->SetRangeUser(0,1);
            Max_Iter_dR[pT-1][cent]->GetYaxis()->SetRangeUser(1e-12,1e-5);

            SetColor(Max_Total_dR[pT-1][cent], 0);
            Max_Total_dR[pT-1][cent]->SetFillStyle(3344);
            Max_Total_dR[pT-1][cent]->SetFillColorAlpha(kGreen-2, 0.6);
            Min_Total_dR[pT-1][cent]->SetFillStyle(1001);
            Min_Total_dR[pT-1][cent]->SetFillColor(10);
            Max_Total_dR[pT-1][cent]->Draw("HIST ][");
            Min_Total_dR[pT-1][cent]->Draw("SAME ][");

            AvgJetdR[pT-1][cent]->SetMarkerStyle(24);
            AvgJetdR[pT-1][cent]->Draw("E1 SAME");
            gPad->RedrawAxis();
        }
    }

    MinMax_dRCanvas->SaveAs(Form("%s/MinMax_dR_Iter_%i.pdf", PlotDir.Data(),4));

    TFile *OldQMdRResults = new TFile("/Volumes/WorkDrive/work/2022/PreliminaryPlots/April2/DeltaR.root", "READ");
    TH1D *DeltaR_Central = (TH1D *)OldQMdRResults->Get("DR_0_10_Sys");
    TH1D *DeltaR_MidCentral = (TH1D *)OldQMdRResults->Get("DR_10_40_Sys");
    TH1D *DeltaR_Peripheral = (TH1D *)OldQMdRResults->Get("DR_40_80_Sys");

    auto QMdRlegend = new TLegend(0.1, 0.1, 0.4, 0.4);
    QMdRlegend->AddEntry(DeltaR_Central, "QM 2022 Version", "lp");
    QMdRlegend->AddEntry(Max_Total_dR[4][0], "Current Version", "f");

    TCanvas *MinMax_dRCanvas_WithOldQMResults = (TCanvas *)MinMax_dRCanvas->DrawClone();
    MinMax_dRCanvas_WithOldQMResults->cd(5);
    SetColor(DeltaR_Central, kRed, 20);
    DeltaR_Central->Draw("EP SAME");
    QMdRlegend->Draw("SAME");
    MinMax_dRCanvas_WithOldQMResults->cd(10);
    SetColor(DeltaR_MidCentral, kRed, 20);
    DeltaR_MidCentral->Draw("EP SAME");
    QMdRlegend->Draw("SAME");
    MinMax_dRCanvas_WithOldQMResults->cd(15);
    SetColor(DeltaR_Peripheral, kRed, 20);
    DeltaR_Peripheral->Draw("EP SAME");
    QMdRlegend->Draw("SAME");

    MinMax_dRCanvas_WithOldQMResults->SaveAs(Form("%s/MinMax_dR_Iter_%i_QM_Comparison.pdf", PlotDir.Data(),4));

    // Now we need to make the total systematic uncertainty for RCPs

    TH1D *RCP_Min_Iter_Pt[5][3];
    TH1D *RCP_Min_Iter_Z[5][3];
    TH1D *RCP_Min_Iter_ZPerJet[5][3];
    TH1D *RCP_Min_Iter_dR[5][3];

    TH1D *RCP_Max_Iter_Pt[5][3];
    TH1D *RCP_Max_Iter_Z[5][3];
    TH1D *RCP_Max_Iter_ZPerJet[5][3];
    TH1D *RCP_Max_Iter_dR[5][3];

    TH1D *RCP_Min_Prior_Pt[5][3];
    TH1D *RCP_Min_Prior_Z[5][3];
    TH1D *RCP_Min_Prior_ZPerJet[5][3];
    TH1D *RCP_Min_Prior_dR[5][3];

    TH1D *RCP_Max_Prior_Pt[5][3];
    TH1D *RCP_Max_Prior_Z[5][3];
    TH1D *RCP_Max_Prior_ZPerJet[5][3];
    TH1D *RCP_Max_Prior_dR[5][3];

    TH1D *RCP_Min_Total_Pt[5][3];
    TH1D *RCP_Min_Total_Z[5][3];
    TH1D *RCP_Min_Total_ZPerJet[5][3];
    TH1D *RCP_Min_Total_dR[5][3];

    TH1D *RCP_Max_Total_Pt[5][3];
    TH1D *RCP_Max_Total_Z[5][3];
    TH1D *RCP_Max_Total_ZPerJet[5][3];
    TH1D *RCP_Max_Total_dR[5][3];

    // for (int pT = 1; pT <= 5; pT++){
    //     for (int cent = 0; cent < 2; cent++){
    //         RCP_Min_Iter_Pt[pT-1][cent] = (TH1D *)Min_Iter_Pt[pT-1][cent]->Clone();
    //         RCP_Max_Iter_Pt[pT-1][cent] = (TH1D *)Max_Iter_Pt[pT-1][cent]->Clone();

    //         SetName(RCP_Min_Iter_Pt[pT-1][cent], Form("%s Min With Iter Sys", RCP_Pt[pT-1][0][cent]->GetName()));
    //         SetName(RCP_Max_Iter_Pt[pT-1][cent], Form("%s Max With Iter Sys", RCP_Pt[pT-1][0][cent]->GetName()));

    //         RCP_Min_Iter_Pt[pT-1][cent]->Divide(Min_Iter_Pt[pT-1][2]);
    //         RCP_Max_Iter_Pt[pT-1][cent]->Divide(Max_Iter_Pt[pT-1][2]);

    //         RCP_Min_Iter_Pt[pT-1][cent]->Scale(taa[2]/taa[cent]);
    //         RCP_Max_Iter_Pt[pT-1][cent]->Scale(taa[2]/taa[cent]);

    //         RCP_Min_Prior_Pt[pT-1][cent] = (TH1D *)Min_Prior_Pt[pT-1][cent]->Clone();
    //         RCP_Max_Prior_Pt[pT-1][cent] = (TH1D *)Max_Prior_Pt[pT-1][cent]->Clone();

    //         SetName(RCP_Min_Prior_Pt[pT-1][cent], Form("%s Min With Prior Sys", RCP_Pt[pT-1][0][cent]->GetName()));
    //         SetName(RCP_Max_Prior_Pt[pT-1][cent], Form("%s Max With Prior Sys", RCP_Pt[pT-1][0][cent]->GetName()));

    //         RCP_Min_Prior_Pt[pT-1][cent]->Divide(Min_Prior_Pt[pT-1][2]);
    //         RCP_Max_Prior_Pt[pT-1][cent]->Divide(Max_Prior_Pt[pT-1][2]);

    //         RCP_Min_Prior_Pt[pT-1][cent]->Scale(taa[2]/taa[cent]);
    //         RCP_Max_Prior_Pt[pT-1][cent]->Scale(taa[2]/taa[cent]);

    //         RCP_Min_Total_Pt[pT-1][cent] = (TH1D *)Min_Total_Pt[pT-1][cent]->Clone();
    //         RCP_Max_Total_Pt[pT-1][cent] = (TH1D *)Max_Total_Pt[pT-1][cent]->Clone();

    //         SetName(RCP_Min_Total_Pt[pT-1][cent], Form("%s Min With Total Sys", RCP_Pt[pT-1][0][cent]->GetName()));
    //         SetName(RCP_Max_Total_Pt[pT-1][cent], Form("%s Max With Total Sys", RCP_Pt[pT-1][0][cent]->GetName()));

    //         RCP_Min_Total_Pt[pT-1][cent]->Divide(Min_Total_Pt[pT-1][2]);
    //         RCP_Max_Total_Pt[pT-1][cent]->Divide(Max_Total_Pt[pT-1][2]);

    //         RCP_Min_Total_Pt[pT-1][cent]->Scale(taa[2]/taa[cent]);
    //         RCP_Max_Total_Pt[pT-1][cent]->Scale(taa[2]/taa[cent]);
    //     }
    // }

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){
            RCP_Min_Iter_Pt[pT-1][cent] = (TH1D *)RCP_Pt[pT-1][0][cent]->Clone();
            RCP_Max_Iter_Pt[pT-1][cent] = (TH1D *)RCP_Pt[pT-1][0][cent]->Clone();
            SetName(RCP_Min_Iter_Pt[pT-1][cent], Form("%s Min With Iter Sys", RCP_Pt[pT-1][0][cent]->GetName()));
            SetName(RCP_Max_Iter_Pt[pT-1][cent], Form("%s Max With Iter Sys", RCP_Pt[pT-1][0][cent]->GetName()));

            RCP_Min_Prior_Pt[pT-1][cent] = (TH1D *)RCP_Pt[pT-1][0][cent]->Clone();
            RCP_Max_Prior_Pt[pT-1][cent] = (TH1D *)RCP_Pt[pT-1][0][cent]->Clone();
            SetName(RCP_Min_Prior_Pt[pT-1][cent], Form("%s Min With Prior Sys", RCP_Pt[pT-1][0][cent]->GetName()));
            SetName(RCP_Max_Prior_Pt[pT-1][cent], Form("%s Max With Prior Sys", RCP_Pt[pT-1][0][cent]->GetName()));

            RCP_Min_Total_Pt[pT-1][cent] = (TH1D *)RCP_Pt[pT-1][0][cent]->Clone();
            RCP_Max_Total_Pt[pT-1][cent] = (TH1D *)RCP_Pt[pT-1][0][cent]->Clone();
            SetName(RCP_Min_Total_Pt[pT-1][cent], Form("%s Min With Total Sys", RCP_Pt[pT-1][0][cent]->GetName()));
            SetName(RCP_Max_Total_Pt[pT-1][cent], Form("%s Max With Total Sys", RCP_Pt[pT-1][0][cent]->GetName()));

            for (int i = 1; i <= RCP_Pt[pT-1][0][cent]->GetNbinsX(); i++){
                
                // double minerriter = TMath::Min(0.01*MinSystematicUncertainty_Pt[pT-1][cent]->GetBinContent(i), 0.01*MinSystematicUncertainty_Pt[pT-1][cent]->GetBinContent(i));
                // double maxerriter = TMath::Max(0.01*MaxSystematicUncertainty_Pt[pT-1][cent]->GetBinContent(i), 0.01*MaxSystematicUncertainty_Pt[pT-1][cent]->GetBinContent(i));

                double minerriter = MinSystematicUncertainty_RCP_Pt[pT-1][cent]->GetBinContent(i)/100.;
                double maxerriter = MaxSystematicUncertainty_RCP_Pt[pT-1][cent]->GetBinContent(i)/100.;

                double minvaliter = RCP_Pt[pT-1][0][cent]->GetBinContent(i)*(1.0 - minerriter);
                double maxvaliter = RCP_Pt[pT-1][0][cent]->GetBinContent(i)*(1.0 + maxerriter);

                // double minerrprior = TMath::Min(0.01*MinSystematicUncertainty_Prior_Pt[pT-1][cent]->GetBinContent(i), 0.01*MinSystematicUncertainty_Prior_Pt[pT-1][cent]->GetBinContent(i));
                // double maxerrprior = TMath::Max(0.01*MaxSystematicUncertainty_Prior_Pt[pT-1][cent]->GetBinContent(i), 0.01*MaxSystematicUncertainty_Prior_Pt[pT-1][cent]->GetBinContent(i));

                double minerrprior = MinSystematicUncertainty_Prior_RCP_Pt[pT-1][cent]->GetBinContent(i)/100.;
                double maxerrprior = MaxSystematicUncertainty_Prior_RCP_Pt[pT-1][cent]->GetBinContent(i)/100.;

                double minvalprior = RCP_Pt[pT-1][0][cent]->GetBinContent(i)*(1.0 - minerrprior);
                double maxvalprior = RCP_Pt[pT-1][0][cent]->GetBinContent(i)*(1.0 + maxerrprior);

                double minvaltotal = RCP_Pt[pT-1][0][cent]->GetBinContent(i)*(1.0 - sqrt(pow(minerriter,2) + pow(minerrprior,2)));
                double maxvaltotal = RCP_Pt[pT-1][0][cent]->GetBinContent(i)*(1.0 + sqrt(pow(maxerriter,2) + pow(maxerrprior,2)));

                RCP_Min_Iter_Pt[pT-1][cent]->SetBinContent(i, minvaliter);
                RCP_Min_Iter_Pt[pT-1][cent]->SetBinError(i, 0);

                RCP_Max_Iter_Pt[pT-1][cent]->SetBinContent(i, maxvaliter);
                RCP_Max_Iter_Pt[pT-1][cent]->SetBinError(i, 0);

                RCP_Min_Prior_Pt[pT-1][cent]->SetBinContent(i, minvalprior);
                RCP_Min_Prior_Pt[pT-1][cent]->SetBinError(i, 0);

                RCP_Max_Prior_Pt[pT-1][cent]->SetBinContent(i, maxvalprior);
                RCP_Max_Prior_Pt[pT-1][cent]->SetBinError(i, 0);

                RCP_Min_Total_Pt[pT-1][cent]->SetBinContent(i, minvaltotal);
                RCP_Min_Total_Pt[pT-1][cent]->SetBinError(i, 0);

                RCP_Max_Total_Pt[pT-1][cent]->SetBinContent(i, maxvaltotal);
                RCP_Max_Total_Pt[pT-1][cent]->SetBinError(i, 0);
            }
        }
    }

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){
            RCP_Min_Iter_Z[pT-1][cent] = (TH1D *)RCP_Z[pT-1][0][cent]->Clone();
            RCP_Max_Iter_Z[pT-1][cent] = (TH1D *)RCP_Z[pT-1][0][cent]->Clone();
            SetName(RCP_Min_Iter_Z[pT-1][cent], Form("%s Min With Iter Sys", RCP_Z[pT-1][0][cent]->GetName()));
            SetName(RCP_Max_Iter_Z[pT-1][cent], Form("%s Max With Iter Sys", RCP_Z[pT-1][0][cent]->GetName()));

            RCP_Min_Prior_Z[pT-1][cent] = (TH1D *)RCP_Z[pT-1][0][cent]->Clone();
            RCP_Max_Prior_Z[pT-1][cent] = (TH1D *)RCP_Z[pT-1][0][cent]->Clone();
            SetName(RCP_Min_Prior_Z[pT-1][cent], Form("%s Min With Prior Sys", RCP_Z[pT-1][0][cent]->GetName()));
            SetName(RCP_Max_Prior_Z[pT-1][cent], Form("%s Max With Prior Sys", RCP_Z[pT-1][0][cent]->GetName()));

            RCP_Min_Total_Z[pT-1][cent] = (TH1D *)RCP_Z[pT-1][0][cent]->Clone();
            RCP_Max_Total_Z[pT-1][cent] = (TH1D *)RCP_Z[pT-1][0][cent]->Clone();
            SetName(RCP_Min_Total_Z[pT-1][cent], Form("%s Min With Total Sys", RCP_Z[pT-1][0][cent]->GetName()));
            SetName(RCP_Max_Total_Z[pT-1][cent], Form("%s Max With Total Sys", RCP_Z[pT-1][0][cent]->GetName()));

            for (int i = 1; i <= RCP_Z[pT-1][0][cent]->GetNbinsX(); i++){
                
                // double minerriter = TMath::Min(0.01*MinSystematicUncertainty_Z[pT-1][cent]->GetBinContent(i), 0.01*MinSystematicUncertainty_Z[pT-1][cent]->GetBinContent(i));
                // double maxerriter = TMath::Max(0.01*MaxSystematicUncertainty_Z[pT-1][cent]->GetBinContent(i), 0.01*MaxSystematicUncertainty_Z[pT-1][cent]->GetBinContent(i));

                double minerriter = MinSystematicUncertainty_RCP_Z[pT-1][cent]->GetBinContent(i)/100.;
                double maxerriter = MaxSystematicUncertainty_RCP_Z[pT-1][cent]->GetBinContent(i)/100.;

                double minvaliter = RCP_Z[pT-1][0][cent]->GetBinContent(i)*(1.0 - minerriter);
                double maxvaliter = RCP_Z[pT-1][0][cent]->GetBinContent(i)*(1.0 + maxerriter);

                // double minerrprior = TMath::Min(0.01*MinSystematicUncertainty_Prior_Z[pT-1][cent]->GetBinContent(i), 0.01*MinSystematicUncertainty_Prior_Z[pT-1][cent]->GetBinContent(i));
                // double maxerrprior = TMath::Max(0.01*MaxSystematicUncertainty_Prior_Z[pT-1][cent]->GetBinContent(i), 0.01*MaxSystematicUncertainty_Prior_Z[pT-1][cent]->GetBinContent(i));

                double minerrprior = MinSystematicUncertainty_Prior_RCP_Z[pT-1][cent]->GetBinContent(i)/100.;
                double maxerrprior = MaxSystematicUncertainty_Prior_RCP_Z[pT-1][cent]->GetBinContent(i)/100.;

                double minvalprior = RCP_Z[pT-1][0][cent]->GetBinContent(i)*(1.0 - minerrprior);
                double maxvalprior = RCP_Z[pT-1][0][cent]->GetBinContent(i)*(1.0 + maxerrprior);

                double minvaltotal = RCP_Z[pT-1][0][cent]->GetBinContent(i)*(1.0 - sqrt(pow(minerriter,2) + pow(minerrprior,2)));
                double maxvaltotal = RCP_Z[pT-1][0][cent]->GetBinContent(i)*(1.0 + sqrt(pow(maxerriter,2) + pow(maxerrprior,2)));

                RCP_Min_Iter_Z[pT-1][cent]->SetBinContent(i, minvaliter);
                RCP_Min_Iter_Z[pT-1][cent]->SetBinError(i, 0);

                RCP_Max_Iter_Z[pT-1][cent]->SetBinContent(i, maxvaliter);
                RCP_Max_Iter_Z[pT-1][cent]->SetBinError(i, 0);

                RCP_Min_Prior_Z[pT-1][cent]->SetBinContent(i, minvalprior);
                RCP_Min_Prior_Z[pT-1][cent]->SetBinError(i, 0);

                RCP_Max_Prior_Z[pT-1][cent]->SetBinContent(i, maxvalprior);
                RCP_Max_Prior_Z[pT-1][cent]->SetBinError(i, 0);

                RCP_Min_Total_Z[pT-1][cent]->SetBinContent(i, minvaltotal);
                RCP_Min_Total_Z[pT-1][cent]->SetBinError(i, 0);

                RCP_Max_Total_Z[pT-1][cent]->SetBinContent(i, maxvaltotal);
                RCP_Max_Total_Z[pT-1][cent]->SetBinError(i, 0);
            }
        }
    }

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){
            RCP_Min_Iter_ZPerJet[pT-1][cent] = (TH1D *)RCP_ZPerJet[pT-1][0][cent]->Clone();
            RCP_Max_Iter_ZPerJet[pT-1][cent] = (TH1D *)RCP_ZPerJet[pT-1][0][cent]->Clone();
            SetName(RCP_Min_Iter_ZPerJet[pT-1][cent], Form("%s Min With Iter Sys", RCP_ZPerJet[pT-1][0][cent]->GetName()));
            SetName(RCP_Max_Iter_ZPerJet[pT-1][cent], Form("%s Max With Iter Sys", RCP_ZPerJet[pT-1][0][cent]->GetName()));

            RCP_Min_Prior_ZPerJet[pT-1][cent] = (TH1D *)RCP_ZPerJet[pT-1][0][cent]->Clone();
            RCP_Max_Prior_ZPerJet[pT-1][cent] = (TH1D *)RCP_ZPerJet[pT-1][0][cent]->Clone();
            SetName(RCP_Min_Prior_ZPerJet[pT-1][cent], Form("%s Min With Prior Sys", RCP_ZPerJet[pT-1][0][cent]->GetName()));
            SetName(RCP_Max_Prior_ZPerJet[pT-1][cent], Form("%s Max With Prior Sys", RCP_ZPerJet[pT-1][0][cent]->GetName()));

            RCP_Min_Total_ZPerJet[pT-1][cent] = (TH1D *)RCP_ZPerJet[pT-1][0][cent]->Clone();
            RCP_Max_Total_ZPerJet[pT-1][cent] = (TH1D *)RCP_ZPerJet[pT-1][0][cent]->Clone();
            SetName(RCP_Min_Total_ZPerJet[pT-1][cent], Form("%s Min With Total Sys", RCP_ZPerJet[pT-1][0][cent]->GetName()));
            SetName(RCP_Max_Total_ZPerJet[pT-1][cent], Form("%s Max With Total Sys", RCP_ZPerJet[pT-1][0][cent]->GetName()));

            for (int i = 1; i <= RCP_ZPerJet[pT-1][0][cent]->GetNbinsX(); i++){
                
                // double minerriter = TMath::Min(0.01*MinSystematicUncertainty_ZPerJet[pT-1][cent]->GetBinContent(i), 0.01*MinSystematicUncertainty_ZPerJet[pT-1][cent]->GetBinContent(i));
                // double maxerriter = TMath::Max(0.01*MaxSystematicUncertainty_ZPerJet[pT-1][cent]->GetBinContent(i), 0.01*MaxSystematicUncertainty_ZPerJet[pT-1][cent]->GetBinContent(i));

                double minerriter = MinSystematicUncertainty_RCP_ZPerJet[pT-1][cent]->GetBinContent(i)/100.;
                double maxerriter = MaxSystematicUncertainty_RCP_ZPerJet[pT-1][cent]->GetBinContent(i)/100.;

                double minvaliter = RCP_ZPerJet[pT-1][0][cent]->GetBinContent(i)*(1.0 - minerriter);
                double maxvaliter = RCP_ZPerJet[pT-1][0][cent]->GetBinContent(i)*(1.0 + maxerriter);

                // double minerrprior = TMath::Min(0.01*MinSystematicUncertainty_Prior_ZPerJet[pT-1][cent]->GetBinContent(i), 0.01*MinSystematicUncertainty_Prior_ZPerJet[pT-1][cent]->GetBinContent(i));
                // double maxerrprior = TMath::Max(0.01*MaxSystematicUncertainty_Prior_ZPerJet[pT-1][cent]->GetBinContent(i), 0.01*MaxSystematicUncertainty_Prior_ZPerJet[pT-1][cent]->GetBinContent(i));

                double minerrprior = MinSystematicUncertainty_Prior_RCP_ZPerJet[pT-1][cent]->GetBinContent(i)/100.;
                double maxerrprior = MaxSystematicUncertainty_Prior_RCP_ZPerJet[pT-1][cent]->GetBinContent(i)/100.;

                double minvalprior = RCP_ZPerJet[pT-1][0][cent]->GetBinContent(i)*(1.0 - minerrprior);
                double maxvalprior = RCP_ZPerJet[pT-1][0][cent]->GetBinContent(i)*(1.0 + maxerrprior);

                double minvaltotal = RCP_ZPerJet[pT-1][0][cent]->GetBinContent(i)*(1.0 - sqrt(pow(minerriter,2) + pow(minerrprior,2)));
                double maxvaltotal = RCP_ZPerJet[pT-1][0][cent]->GetBinContent(i)*(1.0 + sqrt(pow(maxerriter,2) + pow(maxerrprior,2)));

                RCP_Min_Iter_ZPerJet[pT-1][cent]->SetBinContent(i, minvaliter);
                RCP_Min_Iter_ZPerJet[pT-1][cent]->SetBinError(i, 0);

                RCP_Max_Iter_ZPerJet[pT-1][cent]->SetBinContent(i, maxvaliter);
                RCP_Max_Iter_ZPerJet[pT-1][cent]->SetBinError(i, 0);

                RCP_Min_Prior_ZPerJet[pT-1][cent]->SetBinContent(i, minvalprior);
                RCP_Min_Prior_ZPerJet[pT-1][cent]->SetBinError(i, 0);

                RCP_Max_Prior_ZPerJet[pT-1][cent]->SetBinContent(i, maxvalprior);
                RCP_Max_Prior_ZPerJet[pT-1][cent]->SetBinError(i, 0);

                RCP_Min_Total_ZPerJet[pT-1][cent]->SetBinContent(i, minvaltotal);
                RCP_Min_Total_ZPerJet[pT-1][cent]->SetBinError(i, 0);

                RCP_Max_Total_ZPerJet[pT-1][cent]->SetBinContent(i, maxvaltotal);
                RCP_Max_Total_ZPerJet[pT-1][cent]->SetBinError(i, 0);
            }
        }
    }

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){
            RCP_Min_Iter_dR[pT-1][cent] = (TH1D *)RCP_dR[pT-1][0][cent]->Clone();
            RCP_Max_Iter_dR[pT-1][cent] = (TH1D *)RCP_dR[pT-1][0][cent]->Clone();
            SetName(RCP_Min_Iter_dR[pT-1][cent], Form("%s Min With Iter Sys", RCP_dR[pT-1][0][cent]->GetName()));
            SetName(RCP_Max_Iter_dR[pT-1][cent], Form("%s Max With Iter Sys", RCP_dR[pT-1][0][cent]->GetName()));

            RCP_Min_Prior_dR[pT-1][cent] = (TH1D *)RCP_dR[pT-1][0][cent]->Clone();
            RCP_Max_Prior_dR[pT-1][cent] = (TH1D *)RCP_dR[pT-1][0][cent]->Clone();
            SetName(RCP_Min_Prior_dR[pT-1][cent], Form("%s Min With Prior Sys", RCP_dR[pT-1][0][cent]->GetName()));
            SetName(RCP_Max_Prior_dR[pT-1][cent], Form("%s Max With Prior Sys", RCP_dR[pT-1][0][cent]->GetName()));

            RCP_Min_Total_dR[pT-1][cent] = (TH1D *)RCP_dR[pT-1][0][cent]->Clone();
            RCP_Max_Total_dR[pT-1][cent] = (TH1D *)RCP_dR[pT-1][0][cent]->Clone();
            SetName(RCP_Min_Total_dR[pT-1][cent], Form("%s Min With Total Sys", RCP_dR[pT-1][0][cent]->GetName()));
            SetName(RCP_Max_Total_dR[pT-1][cent], Form("%s Max With Total Sys", RCP_dR[pT-1][0][cent]->GetName()));

            for (int i = 1; i <= RCP_dR[pT-1][0][cent]->GetNbinsX(); i++){
                
                // double minerriter = TMath::Min(0.01*MinSystematicUncertainty_dR[pT-1][cent]->GetBinContent(i), 0.01*MinSystematicUncertainty_dR[pT-1][cent]->GetBinContent(i));
                // double maxerriter = TMath::Max(0.01*MaxSystematicUncertainty_dR[pT-1][cent]->GetBinContent(i), 0.01*MaxSystematicUncertainty_dR[pT-1][cent]->GetBinContent(i));

                double minerriter = MinSystematicUncertainty_RCP_dR[pT-1][cent]->GetBinContent(i)/100.;
                double maxerriter = MaxSystematicUncertainty_RCP_dR[pT-1][cent]->GetBinContent(i)/100.;

                double minvaliter = RCP_dR[pT-1][0][cent]->GetBinContent(i)*(1.0 - minerriter);
                double maxvaliter = RCP_dR[pT-1][0][cent]->GetBinContent(i)*(1.0 + maxerriter);

                // double minerrprior = TMath::Min(0.01*MinSystematicUncertainty_Prior_dR[pT-1][cent]->GetBinContent(i), 0.01*MinSystematicUncertainty_Prior_dR[pT-1][cent]->GetBinContent(i));
                // double maxerrprior = TMath::Max(0.01*MaxSystematicUncertainty_Prior_dR[pT-1][cent]->GetBinContent(i), 0.01*MaxSystematicUncertainty_Prior_dR[pT-1][cent]->GetBinContent(i));

                double minerrprior = MinSystematicUncertainty_Prior_RCP_dR[pT-1][cent]->GetBinContent(i)/100.;
                double maxerrprior = MaxSystematicUncertainty_Prior_RCP_dR[pT-1][cent]->GetBinContent(i)/100.;

                double minvalprior = RCP_dR[pT-1][0][cent]->GetBinContent(i)*(1.0 - minerrprior);
                double maxvalprior = RCP_dR[pT-1][0][cent]->GetBinContent(i)*(1.0 + maxerrprior);

                double minvaltotal = RCP_dR[pT-1][0][cent]->GetBinContent(i)*(1.0 - sqrt(pow(minerriter,2) + pow(minerrprior,2)));
                double maxvaltotal = RCP_dR[pT-1][0][cent]->GetBinContent(i)*(1.0 + sqrt(pow(maxerriter,2) + pow(maxerrprior,2)));

                RCP_Min_Iter_dR[pT-1][cent]->SetBinContent(i, minvaliter);
                RCP_Min_Iter_dR[pT-1][cent]->SetBinError(i, 0);

                RCP_Max_Iter_dR[pT-1][cent]->SetBinContent(i, maxvaliter);
                RCP_Max_Iter_dR[pT-1][cent]->SetBinError(i, 0);

                RCP_Min_Prior_dR[pT-1][cent]->SetBinContent(i, minvalprior);
                RCP_Min_Prior_dR[pT-1][cent]->SetBinError(i, 0);

                RCP_Max_Prior_dR[pT-1][cent]->SetBinContent(i, maxvalprior);
                RCP_Max_Prior_dR[pT-1][cent]->SetBinError(i, 0);

                RCP_Min_Total_dR[pT-1][cent]->SetBinContent(i, minvaltotal);
                RCP_Min_Total_dR[pT-1][cent]->SetBinError(i, 0);

                RCP_Max_Total_dR[pT-1][cent]->SetBinContent(i, maxvaltotal);
                RCP_Max_Total_dR[pT-1][cent]->SetBinError(i, 0);
            }
        }
    }

    TLine *RCP_Pt_Line = (TLine *)GetLineAtOne(RCP_Pt[0][0][0], 5, 20);

    TCanvas *RCP_MinMax_PtCanvas = new TCanvas("RCP_MinMax_Iter_PtCanvas", "RCP_MinMax_Iter_PtCanvas", 2100, 700);
    RCP_MinMax_PtCanvas->Divide(5,2);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){
            RCP_MinMax_PtCanvas->cd(cent*5+pT);
            // gPad->SetLogy();
            RCP_Max_Total_Pt[pT-1][cent]->GetXaxis()->SetRangeUser(5,20);
            RCP_Max_Total_Pt[pT-1][cent]->GetYaxis()->SetRangeUser(0, 2);

            SetColor(RCP_Max_Total_Pt[pT-1][cent], 0);
            SetColor(RCP_Min_Total_Pt[pT-1][cent], 0);
            RCP_Max_Total_Pt[pT-1][cent]->SetFillStyle(3344);
            RCP_Max_Total_Pt[pT-1][cent]->SetFillColorAlpha(kGreen-2, 0.6);
            RCP_Min_Total_Pt[pT-1][cent]->SetFillStyle(1001);
            RCP_Min_Total_Pt[pT-1][cent]->SetFillColor(10);
            RCP_Max_Total_Pt[pT-1][cent]->Draw("HIST ][");
            RCP_Min_Total_Pt[pT-1][cent]->Draw("SAME ][");

            AvgRCP_JetpT[pT-1][cent]->SetMarkerStyle(24);
            AvgRCP_JetpT[pT-1][cent]->Draw("E1 SAME");
            RCP_Pt_Line->Draw("SAME");
            gPad->RedrawAxis();
        }
    }

    RCP_MinMax_PtCanvas->SaveAs(Form("%s/RCP_MinMax_Pt_Iter_%i.pdf", PlotDir.Data(),4));

    TFile *OldQMJetPtRCPResults = new TFile("/Volumes/WorkDrive/work/2022/PreliminaryPlots/April2/JetRCP.root", "READ");
    TH1D *JetRCP_Central = (TH1D *)OldQMJetPtRCPResults->Get("JetRCP_0_10_Sys");
    TH1D *JetRCP_MidCentral = (TH1D *)OldQMJetPtRCPResults->Get("JetRCP_10_40_Sys");

    auto QMJetRCPlegend = new TLegend(0.6, 0.6, 0.9, 0.9);
    QMJetRCPlegend->AddEntry(JetRCP_Central, "QM 2022 Version", "lp");
    QMJetRCPlegend->AddEntry(Max_Total_dR[4][0], "Current Version", "f");

    TCanvas *RCP_MinMax_PtCanvas_WithOldQMResults = (TCanvas *)RCP_MinMax_PtCanvas->DrawClone();

    RCP_MinMax_PtCanvas_WithOldQMResults->cd(5);
    SetColor(JetRCP_Central, kRed, 20);
    JetRCP_Central->Draw("EP SAME");
    QMJetRCPlegend->Draw("SAME");
    RCP_MinMax_PtCanvas_WithOldQMResults->cd(10);
    JetRCP_MidCentral->Draw("EP SAME");
    SetColor(JetRCP_MidCentral, kRed, 20);
    QMJetRCPlegend->Draw("SAME");

    RCP_MinMax_PtCanvas_WithOldQMResults->SaveAs(Form("%s/RCP_MinMax_Pt_Iter_%i_QM_Comparison.pdf", PlotDir.Data(),4));

    // Adding the STAR charged particle RCP results here.
    TFile *ChargedParticleFile = new TFile("../ChargedJet_Plots_STAR.root");
    TGraph *ChargedParticleRCP = (TGraph*)ChargedParticleFile->Get("Table 7/Graph1D_y3");

    auto chargedjetlegend = new TLegend(0.1, 0.6, 0.9, 0.9);

    TCanvas *RCP_MinMax_PtCanvas_WithSTAR = (TCanvas *)RCP_MinMax_PtCanvas->DrawClone();
    RCP_MinMax_PtCanvas_WithSTAR->cd(5);
    ChargedParticleRCP->SetMarkerStyle(29);
    ChargedParticleRCP->SetMarkerSize(1.5);
    ChargedParticleRCP->SetMarkerColor(kRed);
    ChargedParticleRCP->SetLineColor(kRed);
    ChargedParticleRCP->Draw("SAME P");
    chargedjetlegend->AddEntry(ChargedParticleRCP, "STAR Charged Jets R = 0.4, p_{T}^{Lead} > 5.0 GeV/#it{c}, 0-10 %/60-80 %", "P");
    chargedjetlegend->AddEntry(RCP_Max_Total_Pt[4][0], "D^{0} Full Jets R = 0.4, p_{T}^{D^{0}} > 5.0 GeV/#it{c}, 0-10 %/40-80 %", "F");
    chargedjetlegend->Draw("SAME");
    gPad->RedrawAxis();
    gPad->Update();
    
    RCP_MinMax_PtCanvas_WithSTAR->SaveAs(Form("%s/RCP_MinMax_Pt_Iter_%i_STAR_ChargedJets_Comparison.pdf", PlotDir.Data(),4));

    TLine *RCP_Z_Line = (TLine *)GetLineAtOne(RCP_Z[0][0][0]);

    TCanvas *RCP_MinMax_ZCanvas = new TCanvas("RCP_MinMax_Iter_ZCanvas", "RCP_MinMax_Iter_ZCanvas", 2100, 700);
    RCP_MinMax_ZCanvas->Divide(5,2);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){
            RCP_MinMax_ZCanvas->cd(cent*5+pT);
            // gPad->SetLogy();
            RCP_Max_Iter_Z[pT-1][cent]->GetXaxis()->SetRangeUser(0, 1);
            RCP_Max_Iter_Z[pT-1][cent]->GetYaxis()->SetRangeUser(0, 2);

            SetColor(RCP_Max_Total_Z[pT-1][cent], 0);
            SetColor(RCP_Min_Total_Z[pT-1][cent], 0);
            RCP_Max_Total_Z[pT-1][cent]->SetFillStyle(3344);
            RCP_Max_Total_Z[pT-1][cent]->SetFillColorAlpha(kGreen-2, 0.6);
            RCP_Min_Total_Z[pT-1][cent]->SetFillStyle(1001);
            RCP_Min_Total_Z[pT-1][cent]->SetFillColor(10);
            RCP_Max_Total_Z[pT-1][cent]->Draw("HIST ][");
            RCP_Min_Total_Z[pT-1][cent]->Draw("SAME ][");

            AvgRCP_JetZ[pT-1][cent]->SetMarkerStyle(24);
            AvgRCP_JetZ[pT-1][cent]->Draw("E1 SAME");
            RCP_Z_Line->Draw("SAME");
            gPad->RedrawAxis();
        }
    }

    RCP_MinMax_ZCanvas->SaveAs(Form("%s/RCP_MinMax_Z_Iter_%i.pdf", PlotDir.Data(),4));

    TLine *RCP_ZPerJet_Line = (TLine *)GetLineAtOne(RCP_ZPerJet[0][0][0]);

    TCanvas *RCP_MinMax_ZPerJetCanvas = new TCanvas("RCP_MinMax_Iter_ZPerJetCanvas", "RCP_MinMax_Iter_ZPerJetCanvas", 2100, 700);
    RCP_MinMax_ZPerJetCanvas->Divide(5,2);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){
            RCP_MinMax_ZPerJetCanvas->cd(cent*5+pT);
            // gPad->SetLogy();
            RCP_Max_Iter_ZPerJet[pT-1][cent]->GetXaxis()->SetRangeUser(0, 1);
            RCP_Max_Iter_ZPerJet[pT-1][cent]->GetYaxis()->SetRangeUser(0, 2);

            SetColor(RCP_Max_Total_ZPerJet[pT-1][cent], 0);
            SetColor(RCP_Min_Total_ZPerJet[pT-1][cent], 0);
            RCP_Max_Total_ZPerJet[pT-1][cent]->SetFillStyle(3344);
            RCP_Max_Total_ZPerJet[pT-1][cent]->SetFillColorAlpha(kGreen-2, 0.6);
            RCP_Min_Total_ZPerJet[pT-1][cent]->SetFillStyle(1001);
            RCP_Min_Total_ZPerJet[pT-1][cent]->SetFillColor(10);
            RCP_Max_Total_ZPerJet[pT-1][cent]->Draw("HIST ][");
            RCP_Min_Total_ZPerJet[pT-1][cent]->Draw("SAME ][");

            AvgRCP_JetZPerJet[pT-1][cent]->SetMarkerStyle(24);
            AvgRCP_JetZPerJet[pT-1][cent]->Draw("E1 SAME");
            RCP_ZPerJet_Line->Draw("SAME");
            gPad->RedrawAxis();
        }
    }

    RCP_MinMax_ZPerJetCanvas->SaveAs(Form("%s/RCP_MinMax_ZPerJet_Iter_%i.pdf", PlotDir.Data(),4));

    TLine *RCP_dR_Line = (TLine *)GetLineAtOne(RCP_dR[0][0][0], 0, 0.2);

    TCanvas *RCP_MinMax_dRCanvas = new TCanvas("RCP_MinMax_dRCanvas", "RCP_MinMax_dRCanvas", 2100, 700);
    RCP_MinMax_dRCanvas->Divide(5,2);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){
            RCP_MinMax_dRCanvas->cd(cent*5+pT);
            // gPad->SetLogy();
            RCP_Max_Iter_dR[pT-1][cent]->GetXaxis()->SetRangeUser(0,0.2);
            RCP_Max_Iter_dR[pT-1][cent]->GetYaxis()->SetRangeUser(0,2);

            SetColor(RCP_Max_Total_dR[pT-1][cent], 0);
            // SetColor(RCP_Min_Total_dR[pT-1][cent], 0);
            RCP_Max_Total_dR[pT-1][cent]->SetFillStyle(3344);
            RCP_Max_Total_dR[pT-1][cent]->SetFillColorAlpha(kGreen-2, 0.6);
            RCP_Min_Total_dR[pT-1][cent]->SetFillStyle(1001);
            RCP_Min_Total_dR[pT-1][cent]->SetFillColor(10);
            RCP_Max_Total_dR[pT-1][cent]->Draw("HIST ][");
            RCP_Min_Total_dR[pT-1][cent]->Draw("SAME ][");

            AvgRCP_JetdR[pT-1][cent]->SetMarkerStyle(24);
            AvgRCP_JetdR[pT-1][cent]->Draw("E1 SAME");
            RCP_dR_Line->Draw("SAME");
            gPad->RedrawAxis();
        }
    }

    RCP_MinMax_dRCanvas->SaveAs(Form("%s/RCP_MinMax_dR_Iter_%i.pdf", PlotDir.Data(),4));

    TFile *OldQMdRRatioResults = new TFile("/Volumes/WorkDrive/Y2020/STAR/STARAnalysis/2021/2022/PreliminaryPlots/April2/DRRatios.root", "READ");
    TH1D *OldQMdR_Cent = (TH1D*)OldQMdRRatioResults->Get("DR_Ratio_0_10");
    TH1D *OldQMdR_Periph = (TH1D*)OldQMdRRatioResults->Get("DR_Ratio_10_40");
    TH1D *OldQMdR_CentSys = (TH1D*)OldQMdRRatioResults->Get("DR_Ratio_0_10_Sys");
    TH1D *OldQMdR_PeriphSys = (TH1D*)OldQMdRRatioResults->Get("DR_Ratio_10_40_Sys");

    for (int i = 1; i <= OldQMdR_Cent->GetNbinsX(); i++){
        OldQMdR_Cent->SetBinError(i, TMath::Sqrt(pow(OldQMdR_CentSys->GetBinError(i), 2) + pow(OldQMdR_Cent->GetBinError(i), 2)));
        OldQMdR_Periph->SetBinError(i, TMath::Sqrt(pow(OldQMdR_PeriphSys->GetBinError(i), 2) + pow(OldQMdR_Periph->GetBinError(i), 2)));
    }

    auto QMdRRatioLegend = new TLegend(0.1, 0.1, 0.4, 0.4);
    QMdRRatioLegend->AddEntry(OldQMdR_CentSys, "QM 2022 Version", "lp");
    QMdRRatioLegend->AddEntry(RCP_Max_Total_dR[4][0], "Current Version", "f");

    TCanvas *RCP_MinMax_dRCanvas_WithOldQMResults = (TCanvas *)RCP_MinMax_dRCanvas->DrawClone();
    RCP_MinMax_dRCanvas_WithOldQMResults->cd(5);
    SetColor(OldQMdR_CentSys, kRed, 20);
    OldQMdR_CentSys->Draw("EP SAME");
    QMdRRatioLegend->Draw("SAME");
    RCP_MinMax_dRCanvas_WithOldQMResults->cd(10);
    SetColor(OldQMdR_PeriphSys, kRed, 20);
    OldQMdR_PeriphSys->Draw("EP SAME");
    QMdRRatioLegend->Draw("SAME");

    RCP_MinMax_dRCanvas_WithOldQMResults->SaveAs(Form("%s/RCP_MinMax_dR_Iter_%i_QM_Comparison.pdf", PlotDir.Data(),4));

    // This here is Multifold

    TFile *MultifoldResults = new TFile("/Volumes/WorkDrive/STAR-Workspace/D0Analysis/Analysis/AuAu/CombinedAnalysisSuite/MultiFold/Aug13/D0Jets4GeV_FromMultifold.root", "READ");
    int iterforcomparison = 4;
    TString Observable[3] = {"MCJetPt", "MCZ", "MCDeltaR"};

    TH1D *Multifold_JetPt[5][3][3];
    TH1D *Multifold_Z[5][3][3];
    TH1D *Multifold_ZPerJet[5][3][3];
    TH1D *Multifold_dR[5][3][3];

    TH1D *Min_Multifold_JetPt[5][3];
    TH1D *Max_Multifold_JetPt[5][3];
    TH1D *Min_Multifold_Z[5][3];
    TH1D *Max_Multifold_Z[5][3];
    TH1D *Min_Multifold_ZPerJet[5][3];
    TH1D *Max_Multifold_ZPerJet[5][3];
    TH1D *Min_Multifold_dR[5][3];
    TH1D *Max_Multifold_dR[5][3];

    TH1D *RCP_Multifold_JetPt[5];
    TH1D *RCP_Multifold_Z[5];
    TH1D *RCP_Multifold_ZPerJet[5];
    TH1D *RCP_Multifold_dR[5];

    for (int iter = 3; iter <= 5; iter++){
        Multifold_JetPt[3][0][iter-3] = (TH1D*)MultifoldResults->Get(Form("Central_MCJetPt_Iter_%i", 2*iter));
        Multifold_JetPt[3][2][iter-3] = (TH1D*)MultifoldResults->Get(Form("Peripheral_MCJetPt_Iter_%i", 2*iter));

        Multifold_Z[3][0][iter-3] = (TH1D*)MultifoldResults->Get(Form("Central_MCZ_Iter_%i", 2*iter));
        SetName(Multifold_Z[3][0][iter-3], Form("%s CrossSection", Multifold_Z[3][0][iter-3]->GetName()));
        Multifold_Z[3][2][iter-3] = (TH1D*)MultifoldResults->Get(Form("Peripheral_MCZ_Iter_%i", 2*iter));
        SetName(Multifold_Z[3][2][iter-3], Form("%s CrossSection", Multifold_Z[3][2][iter-3]->GetName()));

        Multifold_ZPerJet[3][0][iter-3] = (TH1D*)MultifoldResults->Get(Form("Central_MCZ_Iter_%i", 2*iter));
        SetName(Multifold_ZPerJet[3][0][iter-3], Form("%s PerJet", Multifold_ZPerJet[3][0][iter-3]->GetName()));
        Multifold_ZPerJet[3][2][iter-3] = (TH1D*)MultifoldResults->Get(Form("Peripheral_MCZ_Iter_%i", 2*iter));
        SetName(Multifold_ZPerJet[3][2][iter-3], Form("%s PerJet", Multifold_ZPerJet[3][2][iter-3]->GetName()));

        Multifold_dR[3][0][iter-3] = (TH1D*)MultifoldResults->Get(Form("Central_MCDeltaR_Iter_%i", 2*iter));
        Multifold_dR[3][2][iter-3] = (TH1D*)MultifoldResults->Get(Form("Peripheral_MCDeltaR_Iter_%i", 2*iter));
    }

    cout << "Imported all the multifold histograms" << endl;

    for (int cent = 0; cent < 3; cent++){
        if (cent == 1) continue;
        cout << "Cent = " << Multifold_JetPt[3][cent][iterforcomparison-3]->GetNbinsX() << endl;
        Min_Multifold_JetPt[3][cent] = (TH1D *)Multifold_JetPt[3][cent][iterforcomparison-3]->Clone("Min_Multifold_JetPt_Iter=4");
        Max_Multifold_JetPt[3][cent] = (TH1D *)Multifold_JetPt[3][cent][iterforcomparison-3]->Clone("Max_Multifold_JetPt_Iter=4");
        
        for ( int binx = 1; binx <= Multifold_JetPt[3][cent][iterforcomparison-3]->GetNbinsX(); binx++ ){
            Min_Multifold_JetPt[3][cent]->SetBinContent(binx, TMath::Min(Multifold_JetPt[3][cent][0]->GetBinContent(binx), Multifold_JetPt[3][cent][2]->GetBinContent(binx)));
            Max_Multifold_JetPt[3][cent]->SetBinContent(binx, TMath::Max(Multifold_JetPt[3][cent][0]->GetBinContent(binx), Multifold_JetPt[3][cent][2]->GetBinContent(binx)));
        }
    }

    for (int cent = 0; cent < 3; cent++){
        if (cent == 1) continue;
        Min_Multifold_Z[3][cent] = (TH1D *)Multifold_Z[3][cent][iterforcomparison-3]->Clone("Min_Multifold_Z_Iter=4");
        Max_Multifold_Z[3][cent] = (TH1D *)Multifold_Z[3][cent][iterforcomparison-3]->Clone("Max_Multifold_Z_Iter=4");

        for ( int binx = 1; binx <= Multifold_Z[3][cent][iterforcomparison-3]->GetNbinsX(); binx++ ){
            Min_Multifold_Z[3][cent]->SetBinContent(binx, TMath::Min(Multifold_Z[3][cent][0]->GetBinContent(binx), Multifold_Z[3][cent][0]->GetBinContent(binx)));
            Max_Multifold_Z[3][cent]->SetBinContent(binx, TMath::Max(Multifold_Z[3][cent][2]->GetBinContent(binx), Multifold_Z[3][cent][2]->GetBinContent(binx)));
        }
    }

    for (int cent = 0; cent < 3; cent++){
        if (cent == 1) continue;
        Min_Multifold_ZPerJet[3][cent] = (TH1D *)Multifold_ZPerJet[3][cent][iterforcomparison-3]->Clone("Min_Multifold_ZPerJet_Iter=4");
        Max_Multifold_ZPerJet[3][cent] = (TH1D *)Multifold_ZPerJet[3][cent][iterforcomparison-3]->Clone("Max_Multifold_ZPerJet_Iter=4");

        for ( int binx = 1; binx <= Multifold_ZPerJet[3][cent][iterforcomparison-3]->GetNbinsX(); binx++ ){
            Min_Multifold_ZPerJet[3][cent]->SetBinContent(binx, TMath::Min(Multifold_ZPerJet[3][cent][0]->GetBinContent(binx), Multifold_ZPerJet[3][cent][0]->GetBinContent(binx)));
            Max_Multifold_ZPerJet[3][cent]->SetBinContent(binx, TMath::Max(Multifold_ZPerJet[3][cent][2]->GetBinContent(binx), Multifold_ZPerJet[3][cent][2]->GetBinContent(binx)));
        }
    }

    for (int cent = 0; cent < 3; cent++){
        if (cent == 1) continue;
        Min_Multifold_dR[3][cent] = (TH1D *)Multifold_dR[3][cent][iterforcomparison-3]->Clone("Min_Multifold_dR_Iter=4");
        Max_Multifold_dR[3][cent] = (TH1D *)Multifold_dR[3][cent][iterforcomparison-3]->Clone("Max_Multifold_dR_Iter=4");
        for ( int binx = 1; binx <= Multifold_dR[3][cent][iterforcomparison-3]->GetNbinsX(); binx++ ){
            Min_Multifold_dR[3][cent]->SetBinContent(binx, TMath::Min(Multifold_dR[3][cent][0]->GetBinContent(binx), Multifold_dR[3][cent][0]->GetBinContent(binx)));
            Max_Multifold_dR[3][cent]->SetBinContent(binx, TMath::Max(Multifold_dR[3][cent][2]->GetBinContent(binx), Multifold_dR[3][cent][2]->GetBinContent(binx)));
        }
    }
    cout << "Multifold Results" << endl;

    // cout << Multifold_JetPt[3][0]->Integral() << "\t" << Multifold_JetPt[3][2]->Integral() << endl;
    // cout << Multifold_Z[3][0]->Integral() << "\t" << Multifold_Z[3][2]->Integral() << endl;
    // cout << Multifold_dR[3][0]->Integral() << "\t" << Multifold_dR[3][2]->Integral() << endl;

    TFile *MissedRatioFile = new TFile("MultiFoldInputFile_D0Pt_4.root", "READ");
    TH1D *Missed_JetPt_Ratio[3];
    TH1D *Missed_Z_Ratio[3];
    TH1D *Missed_ZPerJet_Ratio[3];
    TH1D *Missed_dR_Ratio[3];

    for (int i = 0; i < 3; i++){
        Missed_JetPt_Ratio[i] = (TH1D*)MissedRatioFile->Get(Form("RatioPtMiss_%i", i));
        Missed_Z_Ratio[i] = (TH1D*)MissedRatioFile->Get(Form("RatioZMiss_%i", i));
        SetName(Missed_Z_Ratio[i], Form("%s CrossSection", Missed_Z_Ratio[i]->GetName()));
        Missed_ZPerJet_Ratio[i] = (TH1D*)MissedRatioFile->Get(Form("RatioZMiss_%i", i));
        SetName(Missed_ZPerJet_Ratio[i], Form("%s PerJet", Missed_ZPerJet_Ratio[i]->GetName()));
        Missed_dR_Ratio[i] = (TH1D*)MissedRatioFile->Get(Form("RatiodRMiss_%i", i));
    }

    cout << "Missed Ratios" << endl;
    // cout << Missed_JetPt_Ratio[0]->Integral() << "\t" << Missed_JetPt_Ratio[2]->Integral() << endl;
    // cout << Missed_Z_Ratio[0]->Integral() << "\t" << Missed_Z_Ratio[2]->Integral() << endl;
    // cout << Missed_dR_Ratio[0]->Integral() << "\t" << Missed_dR_Ratio[2]->Integral() << endl;

    for (int i = 0; i < 3; i++){
        if (i == 1) continue;
        for (int binx = 1; binx <= Multifold_JetPt[3][i][iterforcomparison-3]->GetNbinsX(); binx++){
            int binnum = Missed_JetPt_Ratio[i]->FindBin(Multifold_JetPt[3][i][iterforcomparison-3]->GetBinCenter(binx));
            cout << binnum << "\t" << Missed_JetPt_Ratio[i]->GetBinContent(binnum) << endl;
            Multifold_JetPt[3][i][iterforcomparison-3]->SetBinContent(binx, Multifold_JetPt[3][i][iterforcomparison-3]->GetBinContent(binx)/(1.-Missed_JetPt_Ratio[i]->GetBinContent(binnum)));
            Multifold_JetPt[3][i][iterforcomparison-3]->SetBinError(binx, Multifold_JetPt[3][i][iterforcomparison-3]->GetBinError(binx)/(1.-Missed_JetPt_Ratio[i]->GetBinContent(binnum)));

            Min_Multifold_JetPt[3][i]->SetBinContent(binx, Min_Multifold_JetPt[3][i]->GetBinContent(binx)/(1.-Missed_JetPt_Ratio[i]->GetBinContent(binnum)));
            Min_Multifold_JetPt[3][i]->SetBinError(binx, Min_Multifold_JetPt[3][i]->GetBinError(binx)/(1.-Missed_JetPt_Ratio[i]->GetBinContent(binnum)));
            Max_Multifold_JetPt[3][i]->SetBinContent(binx, Max_Multifold_JetPt[3][i]->GetBinContent(binx)/(1.-Missed_JetPt_Ratio[i]->GetBinContent(binnum)));
            Max_Multifold_JetPt[3][i]->SetBinError(binx, Max_Multifold_JetPt[3][i]->GetBinError(binx)/(1.-Missed_JetPt_Ratio[i]->GetBinContent(binnum)));
        }

        ProcessSpectra(Multifold_JetPt[3][i][iterforcomparison-3]);
        Multifold_JetPt[3][i][iterforcomparison-3]->Scale(1./nevents[i]);
        SetColor(Multifold_JetPt[3][i][iterforcomparison-3], kRed, 34, 1, kRed);
        ProcessSpectra(Min_Multifold_JetPt[3][i]);
        Min_Multifold_JetPt[3][i]->Scale(1./nevents[i]);
        SetColor(Min_Multifold_JetPt[3][i], kRed, 0, 0, kRed);
        ProcessSpectra(Max_Multifold_JetPt[3][i]);
        SetColor(Max_Multifold_JetPt[3][i], kRed, 0, 0, kRed);
        Max_Multifold_JetPt[3][i]->Scale(1./nevents[i]);

        RCP_Multifold_JetPt[3] = (TH1D *)Multifold_JetPt[3][0][iterforcomparison-3]->Clone(Form("RCP_JetPt_Pt%i", 3));
        RCP_Multifold_JetPt[3]->Divide(Multifold_JetPt[3][2][iterforcomparison-3]);
        RCP_Multifold_JetPt[3]->Scale(taa[2]/taa[0]);

        for (int binx = 1; binx <= Multifold_Z[3][i][iterforcomparison-3]->GetNbinsX(); binx++){
            int binnum = Missed_Z_Ratio[i]->FindBin(Multifold_Z[3][i][iterforcomparison-3]->GetBinCenter(binx));
            Multifold_Z[3][i][iterforcomparison-3]->SetBinContent(binx, Multifold_Z[3][i][iterforcomparison-3]->GetBinContent(binx)/(1.-Missed_Z_Ratio[i]->GetBinContent(binnum)));
            Multifold_Z[3][i][iterforcomparison-3]->SetBinError(binx, Multifold_Z[3][i][iterforcomparison-3]->GetBinError(binx)/(1.-Missed_Z_Ratio[i]->GetBinContent(binnum)));

            Min_Multifold_Z[3][i]->SetBinContent(binx, Min_Multifold_Z[3][i]->GetBinContent(binx)/(1.-Missed_Z_Ratio[i]->GetBinContent(binnum)));
            Min_Multifold_Z[3][i]->SetBinError(binx, Min_Multifold_Z[3][i]->GetBinError(binx)/(1.-Missed_Z_Ratio[i]->GetBinContent(binnum)));
            Max_Multifold_Z[3][i]->SetBinContent(binx, Max_Multifold_Z[3][i]->GetBinContent(binx)/(1.-Missed_Z_Ratio[i]->GetBinContent(binnum)));
            Max_Multifold_Z[3][i]->SetBinError(binx, Max_Multifold_Z[3][i]->GetBinError(binx)/(1.-Missed_Z_Ratio[i]->GetBinContent(binnum)));

        }
        ProcessSpectra(Multifold_Z[3][i][iterforcomparison-3]);
        Multifold_Z[3][i][iterforcomparison-3]->Scale(1./nevents[i]);
        SetColor(Multifold_Z[3][i][iterforcomparison-3], kRed, 34, 1, kRed);
        ProcessSpectra(Min_Multifold_Z[3][i]);
        Min_Multifold_Z[3][i]->Scale(1./nevents[i]);
        SetColor(Min_Multifold_Z[3][i], kRed, 0, 0, kRed);
        ProcessSpectra(Max_Multifold_Z[3][i]);
        SetColor(Max_Multifold_Z[3][i], kRed, 0, 0, kRed);
        Max_Multifold_Z[3][i]->Scale(1./nevents[i]);

        RCP_Multifold_Z[3] = (TH1D *)Multifold_Z[3][0][iterforcomparison-3]->Clone(Form("RCP_Z_Pt%i", 3));
        RCP_Multifold_Z[3]->Divide(Multifold_Z[3][2][iterforcomparison-3]);
        RCP_Multifold_Z[3]->Scale(taa[2]/taa[0]);

        for (int binx = 1; binx <= Multifold_ZPerJet[3][i][iterforcomparison-3]->GetNbinsX(); binx++){
            int binnum = Missed_ZPerJet_Ratio[i]->FindBin(Multifold_ZPerJet[3][i][iterforcomparison-3]->GetBinCenter(binx));
            Multifold_ZPerJet[3][i][iterforcomparison-3]->SetBinContent(binx, Multifold_ZPerJet[3][i][iterforcomparison-3]->GetBinContent(binx)/(1.-Missed_ZPerJet_Ratio[i]->GetBinContent(binnum)));
            Multifold_ZPerJet[3][i][iterforcomparison-3]->SetBinError(binx, Multifold_ZPerJet[3][i][iterforcomparison-3]->GetBinError(binx)/(1.-Missed_ZPerJet_Ratio[i]->GetBinContent(binnum)));

            Min_Multifold_ZPerJet[3][i]->SetBinContent(binx, Min_Multifold_ZPerJet[3][i]->GetBinContent(binx)/(1.-Missed_ZPerJet_Ratio[i]->GetBinContent(binnum)));
            Min_Multifold_ZPerJet[3][i]->SetBinError(binx, Min_Multifold_ZPerJet[3][i]->GetBinError(binx)/(1.-Missed_ZPerJet_Ratio[i]->GetBinContent(binnum)));
            Max_Multifold_ZPerJet[3][i]->SetBinContent(binx, Max_Multifold_ZPerJet[3][i]->GetBinContent(binx)/(1.-Missed_ZPerJet_Ratio[i]->GetBinContent(binnum)));
            Max_Multifold_ZPerJet[3][i]->SetBinError(binx, Max_Multifold_ZPerJet[3][i]->GetBinError(binx)/(1.-Missed_ZPerJet_Ratio[i]->GetBinContent(binnum)));

        }
        Multifold_ZPerJet[3][i][iterforcomparison-3]->Scale(1./Multifold_ZPerJet[3][i][iterforcomparison-3]->Integral());
        DivideByBinWidth(Multifold_ZPerJet[3][i][iterforcomparison-3]);
        SetColor(Multifold_Z[3][i][iterforcomparison-3], kRed, 34, 1, kRed);
        Min_Multifold_ZPerJet[3][i]->Scale(1./Min_Multifold_ZPerJet[3][i]->Integral());
        DivideByBinWidth(Min_Multifold_ZPerJet[3][i]);
        SetColor(Min_Multifold_ZPerJet[3][i], kRed, 0, 0, kRed);
        Max_Multifold_ZPerJet[3][i]->Scale(1./Max_Multifold_ZPerJet[3][i]->Integral());
        DivideByBinWidth(Max_Multifold_ZPerJet[3][i]);
        SetColor(Max_Multifold_Z[3][i], kRed, 0, 0, kRed);

        RCP_Multifold_ZPerJet[3] = (TH1D *)Multifold_ZPerJet[3][0][iterforcomparison-3]->Clone(Form("RCP_ZPerJet_Pt%i", 3));
        RCP_Multifold_ZPerJet[3]->Divide(Multifold_ZPerJet[3][2][iterforcomparison-3]);

        for (int binx = 1; binx <= Multifold_dR[3][i][iterforcomparison-3]->GetNbinsX(); binx++){
            int binnum = Missed_dR_Ratio[i]->FindBin(Multifold_dR[3][i][iterforcomparison-3]->GetBinCenter(binx));
            Multifold_dR[3][i][iterforcomparison-3]->SetBinContent(binx, Multifold_dR[3][i][iterforcomparison-3]->GetBinContent(binx)/(1.-Missed_dR_Ratio[i]->GetBinContent(binnum)));
            Multifold_dR[3][i][iterforcomparison-3]->SetBinError(binx, Multifold_dR[3][i][iterforcomparison-3]->GetBinError(binx)/(1.-Missed_dR_Ratio[i]->GetBinContent(binnum)));

            Min_Multifold_dR[3][i]->SetBinContent(binx, Min_Multifold_dR[3][i]->GetBinContent(binx)/(1.-Missed_dR_Ratio[i]->GetBinContent(binnum)));
            Min_Multifold_dR[3][i]->SetBinError(binx, Min_Multifold_dR[3][i]->GetBinError(binx)/(1.-Missed_dR_Ratio[i]->GetBinContent(binnum)));
            Max_Multifold_dR[3][i]->SetBinContent(binx, Max_Multifold_dR[3][i]->GetBinContent(binx)/(1.-Missed_dR_Ratio[i]->GetBinContent(binnum)));
            Max_Multifold_dR[3][i]->SetBinError(binx, Max_Multifold_dR[3][i]->GetBinError(binx)/(1.-Missed_dR_Ratio[i]->GetBinContent(binnum)));
        }
        Multifold_dR[3][i][iterforcomparison-3]->Scale(1./Multifold_dR[3][i][iterforcomparison-3]->Integral());
        DivideByBinWidth(Multifold_dR[3][i][iterforcomparison-3]);
        SetColor(Multifold_dR[3][i][iterforcomparison-3], kRed, 34, 1, kRed);
        Min_Multifold_dR[3][i]->Scale(1./Min_Multifold_dR[3][i]->Integral());
        DivideByBinWidth(Min_Multifold_dR[3][i]);
        SetColor(Min_Multifold_dR[3][i], kRed, 0, 0, kRed);
        Max_Multifold_dR[3][i]->Scale(1./Max_Multifold_dR[3][i]->Integral());
        DivideByBinWidth(Max_Multifold_dR[3][i]);
        SetColor(Max_Multifold_dR[3][i], kRed, 0, 0, kRed);

        RCP_Multifold_dR[3] = (TH1D *)Multifold_dR[3][0][iterforcomparison-3]->Clone(Form("RCP_dR_Pt%i", 3));
        RCP_Multifold_dR[3]->Divide(Multifold_dR[3][2][iterforcomparison-3]);
    }

    cout << "Made the histograms" << endl;

    auto multifoldptlegend = new TLegend(0.1, 0.1, 0.4, 0.4);

    TCanvas *MultifoldCanvas_JetPt = new TCanvas("MultifoldCanvas_JetPt", "MultifoldCanvas_JetPt", 2100, 700);
    MultifoldCanvas_JetPt->Divide(5,3);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            MultifoldCanvas_JetPt->cd(cent*5+pT);
            gPad->SetLogy();
            Max_Iter_Pt[pT-1][cent]->GetXaxis()->SetRangeUser(5,20);
            Max_Iter_Pt[pT-1][cent]->GetYaxis()->SetRangeUser(1e-10,1e-3);

            SetColor(Max_Iter_Pt[pT-1][cent], 0);
            Max_Iter_Pt[pT-1][cent]->SetFillStyle(1001);
            Max_Iter_Pt[pT-1][cent]->SetFillColorAlpha(kGreen-2, 0.6);
            Min_Iter_Pt[pT-1][cent]->SetFillStyle(1001);
            Min_Iter_Pt[pT-1][cent]->SetFillColor(10);
            Max_Iter_Pt[pT-1][cent]->Draw("HIST ][");
            Min_Iter_Pt[pT-1][cent]->Draw("HIST ][ SAME");

            JetpT[pT-1][0][cent]->SetMarkerStyle(24);
            JetpT[pT-1][0][cent]->GetYaxis()->SetRangeUser(1e-10,1e-3);
            JetpT[pT-1][0][cent]->Draw("E1 SAME");

            if (pT == 4 && cent != 1){
                Max_Multifold_JetPt[pT-1][cent]->SetFillStyle(3344);
                // Max_Multifold_JetPt[pT-1][cent]->SetFillColorAlpha(kMagenta, 1.);
                Max_Multifold_JetPt[pT-1][cent]->SetLineWidth(0);
                Max_Multifold_JetPt[pT-1][cent]->SetFillStyle(1001);
                Min_Multifold_JetPt[pT-1][cent]->SetFillColor(10);
                // Min_Multifold_JetPt[pT-1][cent]->SetLineWidth(0);
                Max_Multifold_JetPt[pT-1][cent]->Draw("HIST ][ SAME");
                Min_Multifold_JetPt[pT-1][cent]->Draw("HIST SAME ][");

                Multifold_JetPt[pT-1][cent][iterforcomparison-3]->Draw("E5 SAME");
                if (cent == 0){
                    multifoldptlegend->AddEntry(JetpT[pT-1][0][cent], "Binned IBU Iter=4", "lep");
                    multifoldptlegend->AddEntry(Multifold_JetPt[pT-1][cent][iterforcomparison-3], "Multifold Iter=4", "lep");
                }
                multifoldptlegend->Draw("SAME");
            }

            gPad->RedrawAxis();
        }
    }

    MultifoldCanvas_JetPt->SaveAs(Form("%s/MultifoldCanvas_JetPt_Iter_%i.pdf", PlotDir.Data(),4));

    auto multifoldzlegend = new TLegend(0.1, 0.6, 0.4, 0.9);

    TCanvas *MultifoldCanvas_Z = new TCanvas("MultifoldCanvas_Z", "MultifoldCanvas_Z", 2100, 700);
    MultifoldCanvas_Z->Divide(5,3);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            MultifoldCanvas_Z->cd(cent*5+pT);
            gPad->SetLogy();
            Max_Iter_Z[pT-1][cent]->GetXaxis()->SetRangeUser(0,0.5);
            Max_Iter_Z[pT-1][cent]->GetYaxis()->SetRangeUser(1e-8,1e-1);

            SetColor(Max_Iter_Z[pT-1][cent], 0);
            Max_Iter_Z[pT-1][cent]->SetFillStyle(1001);
            Max_Iter_Z[pT-1][cent]->SetFillColorAlpha(kGreen-2, 0.6);
            Min_Iter_Z[pT-1][cent]->SetFillStyle(1001);
            Min_Iter_Z[pT-1][cent]->SetFillColor(10);
            Max_Iter_Z[pT-1][cent]->Draw("HIST ][");
            Min_Iter_Z[pT-1][cent]->Draw("HIST ][ SAME");

            JetZ[pT-1][0][cent]->SetMarkerStyle(24);
            JetZ[pT-1][0][cent]->GetYaxis()->SetRangeUser(1e-8,1e-1);
            JetZ[pT-1][0][cent]->Draw("E1");

            if (pT == 4 && cent != 1){
                Max_Multifold_Z[pT-1][cent]->SetFillStyle(3344);
                Max_Multifold_Z[pT-1][cent]->SetLineWidth(0);
                // Max_Multifold_Z[pT-1][cent]->SetFillColorAlpha(kRed, 0.5);
                Min_Multifold_Z[pT-1][cent]->SetFillStyle(1001);
                Min_Multifold_Z[pT-1][cent]->SetFillColor(10);
                Max_Multifold_Z[pT-1][cent]->Draw("HIST ][ SAME");
                Min_Multifold_Z[pT-1][cent]->Draw("HIST SAME ][");

                Multifold_Z[pT-1][cent][iterforcomparison-3]->Draw("E5 SAME");
                if (cent == 0){
                    multifoldzlegend->AddEntry(JetZ[pT-1][0][cent], "Binned IBU Iter=4", "lep");
                    multifoldzlegend->AddEntry(Multifold_Z[pT-1][cent][iterforcomparison-3], "Multifold Iter=4", "lep");
                }
                multifoldzlegend->Draw("SAME");
            }

            gPad->RedrawAxis();
        }
    }

    MultifoldCanvas_Z->SaveAs(Form("%s/MultifoldCanvas_Z_Iter_%i.pdf", PlotDir.Data(),4));

    auto multifoldzperjetlegend = new TLegend(0.1, 0.6, 0.4, 0.9);

    TCanvas *MultifoldCanvas_ZPerJet = new TCanvas("MultifoldCanvas_ZPerJet", "MultifoldCanvas_ZPerJet", 2100, 700);
    MultifoldCanvas_ZPerJet->Divide(5,3);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            MultifoldCanvas_ZPerJet->cd(cent*5+pT);
            gPad->SetLogy();
            Max_Iter_ZPerJet[pT-1][cent]->GetXaxis()->SetRangeUser(0,0.5);
            Max_Iter_ZPerJet[pT-1][cent]->GetYaxis()->SetRangeUser(1e-4,1e4);

            SetColor(Max_Iter_ZPerJet[pT-1][cent], 0);
            Max_Iter_ZPerJet[pT-1][cent]->SetFillStyle(1001);
            Max_Iter_ZPerJet[pT-1][cent]->SetFillColorAlpha(kGreen-2, 0.6);
            Min_Iter_ZPerJet[pT-1][cent]->SetFillStyle(1001);
            Min_Iter_ZPerJet[pT-1][cent]->SetFillColor(10);
            Max_Iter_ZPerJet[pT-1][cent]->Draw("HIST ][");
            Min_Iter_ZPerJet[pT-1][cent]->Draw("HIST ][ SAME");

            JetZPerJet[pT-1][0][cent]->SetMarkerStyle(24);
            JetZPerJet[pT-1][0][cent]->GetYaxis()->SetRangeUser(1e-4,1e4);
            JetZPerJet[pT-1][0][cent]->Draw("E1");

            if (pT == 4 && cent != 1){
                Max_Multifold_ZPerJet[pT-1][cent]->SetFillStyle(3344);
                Max_Multifold_ZPerJet[pT-1][cent]->SetLineWidth(0);
                // Max_Multifold_Z[pT-1][cent]->SetFillColorAlpha(kRed, 0.5);
                Min_Multifold_ZPerJet[pT-1][cent]->SetFillStyle(1001);
                Min_Multifold_ZPerJet[pT-1][cent]->SetFillColor(10);
                Max_Multifold_ZPerJet[pT-1][cent]->Draw("HIST ][ SAME");
                Min_Multifold_ZPerJet[pT-1][cent]->Draw("HIST SAME ][");

                Multifold_ZPerJet[pT-1][cent][iterforcomparison-3]->Draw("E5 SAME");
                if (cent == 0){
                    multifoldzperjetlegend->AddEntry(JetZPerJet[pT-1][0][cent], "Binned IBU Iter=4", "lep");
                    multifoldzperjetlegend->AddEntry(Multifold_ZPerJet[pT-1][cent][iterforcomparison-3], "Multifold Iter=4", "lep");
                }
                multifoldzperjetlegend->Draw("SAME");
            }

            gPad->RedrawAxis();
        }
    }

    MultifoldCanvas_ZPerJet->SaveAs(Form("%s/MultifoldCanvas_ZPerJet_Iter_%i.pdf", PlotDir.Data(),4));

    auto multifolddrlegend = new TLegend(0.6, 0.6, 0.9, 0.9);

    TCanvas *MultifoldCanvas_dR = new TCanvas("MultifoldCanvas_dR", "MultifoldCanvas_dR", 2100, 700);
    MultifoldCanvas_dR->Divide(5,3);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            MultifoldCanvas_dR->cd(cent*5+pT);
            gPad->SetLogy();
            Max_Iter_dR[pT-1][cent]->GetXaxis()->SetRangeUser(0,0.5);
            Max_Iter_dR[pT-1][cent]->GetYaxis()->SetRangeUser(1e-4,1e+4);

            SetColor(Max_Iter_dR[pT-1][cent], 0);
            Max_Iter_dR[pT-1][cent]->SetFillStyle(1001);
            Max_Iter_dR[pT-1][cent]->SetLineWidth(0);
            Max_Iter_dR[pT-1][cent]->SetFillColorAlpha(kGreen-2, 0.6);
            Min_Iter_dR[pT-1][cent]->SetFillStyle(1001);
            Min_Iter_dR[pT-1][cent]->SetFillColor(10);
            Max_Iter_dR[pT-1][cent]->Draw("HIST ][");
            Min_Iter_dR[pT-1][cent]->Draw("HIST ][ SAME");
            
            JetdR[pT-1][0][cent]->SetMarkerStyle(24);
            JetdR[pT-1][0][cent]->GetYaxis()->SetRangeUser(1e-4,1e+4);
            JetdR[pT-1][0][cent]->Draw("E1");

            if (pT == 4 && cent != 1){
                Max_Multifold_dR[pT-1][cent]->SetFillStyle(3344);
                Max_Multifold_dR[pT-1][cent]->SetFillColorAlpha(kRed, 0.5);
                Max_Multifold_dR[pT-1][cent]->SetFillStyle(1001);
                Min_Multifold_dR[pT-1][cent]->SetFillColor(10);
                Max_Multifold_dR[pT-1][cent]->Draw("][ SAME");
                Min_Multifold_dR[pT-1][cent]->Draw("SAME ][");

                Multifold_dR[pT-1][cent][iterforcomparison-3]->Draw("E5 SAME");
                if (cent == 0){
                    multifolddrlegend->AddEntry(JetdR[pT-1][0][cent], "Binned IBU Iter=4", "lep");
                    multifolddrlegend->AddEntry(Multifold_dR[pT-1][cent][iterforcomparison-3], "Multifold Iter=4", "lep");
                }
                multifolddrlegend->Draw("SAME");
            }

            gPad->RedrawAxis();
        }
    }

    MultifoldCanvas_dR->SaveAs(Form("%s/MultifoldCanvas_dR_Iter_%i.pdf", PlotDir.Data(),4));

    TCanvas *RCP_Multifold_PtCanvas = new TCanvas("RCP_Multifold_PtCanvas", "RCP_Multifold_PtCanvas", 2100, 700);
    RCP_Multifold_PtCanvas->Divide(5,2);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){
            RCP_Multifold_PtCanvas->cd(cent*5+pT);
            // gPad->SetLogy();
            RCP_Max_Iter_Pt[pT-1][cent]->GetXaxis()->SetRangeUser(5,20);
            RCP_Max_Iter_Pt[pT-1][cent]->GetYaxis()->SetRangeUser(0, 2);

            SetColor(RCP_Max_Total_Pt[pT-1][cent], 0);
            SetColor(RCP_Min_Total_Pt[pT-1][cent], 0);
            RCP_Max_Total_Pt[pT-1][cent]->SetFillStyle(3344);
            RCP_Max_Total_Pt[pT-1][cent]->SetFillColorAlpha(kGreen-2, 0.6);
            RCP_Min_Total_Pt[pT-1][cent]->SetFillStyle(1001);
            RCP_Min_Total_Pt[pT-1][cent]->SetFillColor(10);
            RCP_Max_Total_Pt[pT-1][cent]->Draw("HIST ][");
            RCP_Min_Total_Pt[pT-1][cent]->Draw("SAME ][");

            RCP_Pt[pT-1][0][cent]->SetMarkerStyle(24);
            RCP_Pt[pT-1][0][cent]->Draw("E1SAME");

            if (pT == 4 && cent == 0) RCP_Multifold_JetPt[pT-1]->Draw("E5 SAME");

            gPad->RedrawAxis();
        }
    }

    RCP_Multifold_PtCanvas->SaveAs(Form("%s/RCP_Multifold_PtCanvas_Iter_%i.pdf", PlotDir.Data(),4));

    TCanvas *RCP_Multifold_ZCanvas = new TCanvas("RCP_Multifold_ZCanvas", "RCP_Multifold_ZCanvas", 2100, 700);
    RCP_Multifold_ZCanvas->Divide(5,2);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){
            RCP_Multifold_ZCanvas->cd(cent*5+pT);
            // gPad->SetLogy();
            RCP_Max_Iter_Z[pT-1][cent]->GetXaxis()->SetRangeUser(0,0.5);
            RCP_Max_Iter_Z[pT-1][cent]->GetYaxis()->SetRangeUser(0, 2);

            SetColor(RCP_Max_Total_Z[pT-1][cent], 0);
            SetColor(RCP_Min_Total_Z[pT-1][cent], 0);
            RCP_Max_Total_Z[pT-1][cent]->SetFillStyle(3344);
            RCP_Max_Total_Z[pT-1][cent]->SetFillColorAlpha(kGreen-2, 0.6);
            RCP_Min_Total_Z[pT-1][cent]->SetFillStyle(1001);
            RCP_Min_Total_Z[pT-1][cent]->SetFillColor(10);
            RCP_Max_Total_Z[pT-1][cent]->Draw("HIST ][");
            RCP_Min_Total_Z[pT-1][cent]->Draw("SAME ][");

            RCP_Z[pT-1][0][cent]->SetMarkerStyle(24);
            RCP_Z[pT-1][0][cent]->Draw("E1SAME");

            if (pT == 4 && cent == 0) RCP_Multifold_Z[pT-1]->Draw("E5 SAME");

            gPad->RedrawAxis();
        }
    }

    RCP_Multifold_ZCanvas->SaveAs(Form("%s/RCP_Multifold_ZCanvas_Iter_%i.pdf", PlotDir.Data(),4));

    TCanvas *RCP_Multifold_ZPerJetCanvas = new TCanvas("RCP_Multifold_ZPerJetCanvas", "RCP_Multifold_ZPerJetCanvas", 2100, 700);
    RCP_Multifold_ZPerJetCanvas->Divide(5,2);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){
            RCP_Multifold_ZPerJetCanvas->cd(cent*5+pT);
            // gPad->SetLogy();
            RCP_Max_Iter_ZPerJet[pT-1][cent]->GetXaxis()->SetRangeUser(0,0.5);
            RCP_Max_Iter_ZPerJet[pT-1][cent]->GetYaxis()->SetRangeUser(0, 2);

            SetColor(RCP_Max_Total_ZPerJet[pT-1][cent], 0);
            SetColor(RCP_Min_Total_ZPerJet[pT-1][cent], 0);
            RCP_Max_Total_ZPerJet[pT-1][cent]->SetFillStyle(3344);
            RCP_Max_Total_ZPerJet[pT-1][cent]->SetFillColorAlpha(kGreen-2, 0.6);
            RCP_Min_Total_ZPerJet[pT-1][cent]->SetFillStyle(1001);
            RCP_Min_Total_ZPerJet[pT-1][cent]->SetFillColor(10);
            RCP_Max_Total_ZPerJet[pT-1][cent]->Draw("HIST ][");
            RCP_Min_Total_ZPerJet[pT-1][cent]->Draw("SAME ][");

            RCP_ZPerJet[pT-1][0][cent]->SetMarkerStyle(24);
            RCP_ZPerJet[pT-1][0][cent]->Draw("E1SAME");

            if (pT == 4 && cent == 0) RCP_Multifold_ZPerJet[pT-1]->Draw("E5 SAME");

            gPad->RedrawAxis();
        }
    }

    RCP_Multifold_ZPerJetCanvas->SaveAs(Form("%s/RCP_Multifold_ZPerJetCanvas_Iter_%i.pdf", PlotDir.Data(),4));

    TCanvas *RCP_Multifold_dRCanvas = new TCanvas("RCP_Multifold_dRCanvas", "RCP_Multifold_dRCanvas", 2100, 700);
    RCP_Multifold_dRCanvas->Divide(5,2);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){
            RCP_Multifold_dRCanvas->cd(cent*5+pT);
            // gPad->SetLogy();
            RCP_Max_Iter_dR[pT-1][cent]->GetXaxis()->SetRangeUser(0,0.5);
            RCP_Max_Iter_dR[pT-1][cent]->GetYaxis()->SetRangeUser(0, 2);

            SetColor(RCP_Max_Total_dR[pT-1][cent], 0);
            SetColor(RCP_Min_Total_dR[pT-1][cent], 0);
            RCP_Max_Total_dR[pT-1][cent]->SetFillStyle(3344);
            RCP_Max_Total_dR[pT-1][cent]->SetFillColorAlpha(kGreen-2, 0.6);
            RCP_Min_Total_dR[pT-1][cent]->SetFillStyle(1001);
            RCP_Min_Total_dR[pT-1][cent]->SetFillColor(10);
            RCP_Max_Total_dR[pT-1][cent]->Draw("HIST ][");
            RCP_Min_Total_dR[pT-1][cent]->Draw("SAME ][");

            RCP_dR[pT-1][0][cent]->SetMarkerStyle(24);
            RCP_dR[pT-1][0][cent]->Draw("E1SAME");

            if (pT == 4 && cent == 0) RCP_Multifold_dR[pT-1]->Draw("E5 SAME");

            gPad->RedrawAxis();
        }
    }

    // End of Multifold

    // START THEORY HERE

    // Central D0 Jet pT Spectra
    
    TH1D *TheoryD0JetPt[3][5];
    TH1D *TheoryD0JetZ[3][5];
    TH1D *TheoryD0JetZPerJet[3][5];
    TH1D *TheoryD0JetdR[3][5];

    TH1D *TheoryRCP_JetPt[2][5];

    int nz_bins_theory = 9;
    double nbinsz_theory[10] = {0., 0.1, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            TheoryD0JetPt[cent][pT-1] = new TH1D(Form("TheoryD0JetPt_%i_%i", cent, pT), Form("TheoryD0JetPt_%i_%i", cent, pT), nbins_jpt, binning_jpt);
            TheoryD0JetZ[cent][pT-1] = new TH1D(Form("TheoryD0JetZ_%i_%i", cent, pT), Form("TheoryD0JetZ_%i_%i", cent, pT), nz_bins_theory, nbinsz_theory);
            TheoryD0JetZPerJet[cent][pT-1] = new TH1D(Form("TheoryD0JetZPerJet_%i_%i", cent, pT), Form("TheoryD0JetZPerJet_%i_%i", cent, pT), nz_bins_theory, nbinsz_theory);
            TheoryD0JetdR[cent][pT-1] = new TH1D(Form("TheoryD0JetdR_%i_%i", cent, pT), Form("TheoryD0JetdR_%i_%i", cent, pT), ndrbins, drbins);
        }

        for (int cent = 0; cent < 2; cent++){
            TheoryRCP_JetPt[cent][pT-1] = new TH1D(Form("TheoryRCP_JetPt_%i_%i", cent, pT), Form("TheoryRCP_JetPt_%i_%i", cent, pT), nbins_jpt, binning_jpt);
        }
    }

    TH1D *DataD0JetPt[3][5];
    TH1D *DataD0JetZ[3][5];
    TH1D *DataD0JetZPerJet[3][5];
    TH1D *DataD0JetdR[3][5];

    double inelcrosssection = 43.82; // mb

    ReadTheoryFiles("LIDO-projection/0-10/D0-jet-spectra-pTD-1-10.dat", TheoryD0JetPt[0][0], 1./inelcrosssection);
    ReadTheoryFiles("LIDO-projection/0-10/D0-jet-spectra-pTD-4-10.dat", TheoryD0JetPt[0][3], 1./inelcrosssection);
    ReadTheoryFiles("LIDO-projection-AuAu200-0-10/D0-jet-spectra-pTD-5-10.dat", TheoryD0JetPt[0][4], 1./inelcrosssection);

    ReadTheoryFiles("LIDO-projection/40-80/D0-jet-spectra-pTD-1-10.dat", TheoryD0JetPt[2][0], 1./inelcrosssection);
    ReadTheoryFiles("LIDO-projection/40-80/D0-jet-spectra-pTD-4-10.dat", TheoryD0JetPt[2][3], 1./inelcrosssection);

    ReadTheoryFiles("LIDO-projection/0-10/D0inJet_z_cross-section-pTD-1-10.dat", TheoryD0JetZ[0][0], 1./inelcrosssection);
    ReadTheoryFiles("LIDO-projection/0-10/D0inJet_z_cross-section-pTD-4-10.dat", TheoryD0JetZ[0][3], 1./inelcrosssection);
    ReadTheoryFiles("LIDO-projection-AuAu200-0-10/D0inJet_z_cross-section-pTD-5-10.dat", TheoryD0JetZ[0][4], 1./inelcrosssection);

    ReadTheoryFiles("LIDO-projection/40-80/D0inJet_z_cross-section-pTD-1-10.dat", TheoryD0JetZ[2][0], 1./inelcrosssection);
    ReadTheoryFiles("LIDO-projection/40-80/D0inJet_z_cross-section-pTD-4-10.dat", TheoryD0JetZ[2][3], 1./inelcrosssection);

    ReadTheoryFiles("LIDO-projection/0-10/D0inJet_z_per-jet-pTD-1-10.dat", TheoryD0JetZPerJet[0][0], 1.);
    ReadTheoryFiles("LIDO-projection/0-10/D0inJet_z_per-jet-pTD-4-10.dat", TheoryD0JetZPerJet[0][3], 1.);

    ReadTheoryFiles("LIDO-projection/40-80/D0inJet_z_per-jet-pTD-1-10.dat", TheoryD0JetZPerJet[2][0], 1.);
    ReadTheoryFiles("LIDO-projection/40-80/D0inJet_z_per-jet-pTD-4-10.dat", TheoryD0JetZPerJet[2][3], 1.);

    ReadTheoryFiles("LIDO-projection/0-10/D0inJet_r_per-jet-pTD-1-10.dat", TheoryD0JetdR[0][0], 1.);
    ReadTheoryFiles("LIDO-projection/0-10/D0inJet_r_per-jet-pTD-4-10.dat", TheoryD0JetdR[0][3], 1.);

    ReadTheoryFiles("LIDO-projection/40-80/D0inJet_r_per-jet-pTD-1-10.dat", TheoryD0JetdR[2][0], 1.);
    ReadTheoryFiles("LIDO-projection/40-80/D0inJet_r_per-jet-pTD-4-10.dat", TheoryD0JetdR[2][3], 1.);

    ReadRCPTheoryFiles("LIDO-projection/RCP-1-10.dat", TheoryRCP_JetPt[0][0], 1.);
    ReadRCPTheoryFiles("LIDO-projection/RCP-4-10.dat", TheoryRCP_JetPt[0][3], 1.);

    TCanvas *TheoryD0JetPtCanvas = new TCanvas("TheoryD0JetPtCanvas", "TheoryD0JetPtCanvas", 2100, 700);
    TheoryD0JetPtCanvas->Divide(5, 3);

    TLegend *theoryvdatalegend[5][3]; 

    TH1D *Max_Total_Pt_Scaled[5][3];
    TH1D *Min_Total_Pt_Scaled[5][3];

    TH1D *tmp = (TH1D *)TheoryD0JetPt[0][4]->Clone();

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            TheoryD0JetPtCanvas->cd(cent*5+pT);
            theoryvdatalegend[pT-1][cent] = new TLegend(0.4, 0.6, 0.9, 0.9);
            gPad->SetLogy();
            DataD0JetPt[cent][pT-1] = (TH1D*)AvgJetpT[pT-1][cent]->Clone(Form("DataD0JetPt_%i_%i", cent, pT));
            SetName(DataD0JetPt[cent][pT-1], Form("DataD0JetPt_%i_%i", cent, pT));
            DataD0JetPt[cent][pT-1]->Scale(2*TMath::Pi()/taa[cent]);
            Max_Total_Pt_Scaled[pT-1][cent] = (TH1D*)Max_Total_Pt[pT-1][cent]->Clone(Form("Max_Total_Pt_Scaled_%i_%i", cent, pT));
            SetName(Max_Total_Pt_Scaled[pT-1][cent], Form("Max_Total_Pt_Scaled_%i_%i", cent, pT));
            Max_Total_Pt_Scaled[pT-1][cent]->Scale(2*TMath::Pi()/taa[cent]);
            Min_Total_Pt_Scaled[pT-1][cent] = (TH1D*)Min_Total_Pt[pT-1][cent]->Clone(Form("Min_Total_Pt_Scaled_%i_%i", cent, pT));
            SetName(Min_Total_Pt_Scaled[pT-1][cent], Form("Min_Total_Pt_Scaled_%i_%i", cent, pT));
            Min_Total_Pt_Scaled[pT-1][cent]->Scale(2*TMath::Pi()/taa[cent]);
            Max_Total_Pt_Scaled[pT-1][cent]->GetYaxis()->SetRangeUser(1e-12, 1e-4);
            Max_Total_Pt_Scaled[pT-1][cent]->Draw("HIST ][");
            Min_Total_Pt_Scaled[pT-1][cent]->Draw("SAME ][");
            DataD0JetPt[cent][pT-1]->Draw("E1 SAME");
            SetColor(TheoryD0JetPt[cent][pT-1], kRed, 20, 1, kRed);
            TheoryD0JetPt[cent][pT-1]->Draw("E5 SAME");
            // SetColor(tmp, kBlack, 20, 1, kBlack);
            // tmp->Draw("E5 SAME");
            theoryvdatalegend[pT-1][cent]->SetHeader(Form("%i < D0 pT [GeV/#it{c}] < %i", pT, 10), "C");
            theoryvdatalegend[pT-1][cent]->AddEntry(DataD0JetPt[cent][pT-1], "STAR AuAu 0-10%", "P");
            theoryvdatalegend[pT-1][cent]->AddEntry(TheoryD0JetPt[cent][pT-1], "LIDO AuAu 0-10% (MPI = off)", "LP");
            theoryvdatalegend[pT-1][cent]->Draw("SAME");
            gPad->RedrawAxis();
        }
    }

    TheoryD0JetPtCanvas->SaveAs(Form("%s/TheoryD0JetPtCanvas.pdf", PlotDir.Data()));

    TLegend *theoryvdatazlegend[5][3];

    TH1D *Max_Total_Z_Scaled[5][3];
    TH1D *Min_Total_Z_Scaled[5][3];

    TCanvas *TheoryD0JetZCanvas = new TCanvas("TheoryD0JetZCanvas", "TheoryD0JetZCanvas", 2100, 700);
    TheoryD0JetZCanvas->Divide(5, 3);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            TheoryD0JetZCanvas->cd(cent*5+pT);
            theoryvdatazlegend[pT-1][cent] = new TLegend(0.1, 0.6, 0.5, 0.9);
            gPad->SetLogy();
            DataD0JetZ[cent][pT-1] = (TH1D*)AvgJetZ[pT-1][cent]->Clone(Form("DataD0JetZ_%i_%i", cent, pT));
            SetName(DataD0JetZ[cent][pT-1], Form("DataD0JetZ_%i_%i", cent, pT));
            DataD0JetZ[cent][pT-1]->Scale(2*TMath::Pi()/taa[cent]);
            Max_Total_Z_Scaled[pT-1][cent] = (TH1D*)Max_Total_Z[pT-1][cent]->Clone(Form("Max_Total_Z_Scaled_%i_%i", cent, pT));
            SetName(Max_Total_Z_Scaled[pT-1][cent], Form("Max_Total_Z_Scaled_%i_%i", cent, pT));
            Max_Total_Z_Scaled[pT-1][cent]->Scale(2*TMath::Pi()/taa[cent]);
            Min_Total_Z_Scaled[pT-1][cent] = (TH1D*)Min_Total_Z[pT-1][cent]->Clone(Form("Min_Total_Z_Scaled_%i_%i", cent, pT));
            SetName(Min_Total_Z_Scaled[pT-1][cent], Form("Min_Total_Z_Scaled_%i_%i", cent, pT));
            Max_Total_Z_Scaled[pT-1][cent]->GetYaxis()->SetRangeUser(1e-8, 1e-1);
            Min_Total_Z_Scaled[pT-1][cent]->Scale(2*TMath::Pi()/taa[cent]);
            Max_Total_Z_Scaled[pT-1][cent]->Draw("HIST ][");
            Min_Total_Z_Scaled[pT-1][cent]->Draw("SAME ][");
            DataD0JetZ[cent][pT-1]->Draw("E1 SAME");
            SetColor(TheoryD0JetZ[cent][pT-1], kRed, 20, 1, kRed);
            TheoryD0JetZ[cent][pT-1]->Draw("E5 SAME");
            theoryvdatazlegend[pT-1][cent]->SetHeader(Form("%i < D0 pT [GeV/#it{c}] < %i", pT, 10), "C");
            theoryvdatazlegend[pT-1][cent]->AddEntry(DataD0JetZ[cent][pT-1], "STAR AuAu 0-10%", "P");
            theoryvdatazlegend[pT-1][cent]->AddEntry(TheoryD0JetZ[cent][pT-1], "LIDO AuAu 0-10% (MPI = off)", "LP");
            theoryvdatazlegend[pT-1][cent]->Draw("SAME");
            gPad->RedrawAxis();
        }
    }

    TheoryD0JetZCanvas->SaveAs(Form("%s/TheoryD0JetZCanvas.pdf", PlotDir.Data()));

    // TLegend *theoryvdatazperjetlegend[5][3];

    // TCanvas *TheoryD0JetZPerJetCanvas = new TCanvas("TheoryD0JetZPerJetCanvas", "TheoryD0JetZPerJetCanvas", 2100, 700);
    // TheoryD0JetZPerJetCanvas->Divide(5, 3);

    // for (int pT = 1; pT <= 5; pT++){
    //     for (int cent = 0; cent < 3; cent++){
    //         TheoryD0JetZPerJetCanvas->cd(cent*5+pT);
    //         theoryvdatazperjetlegend[pT-1][cent] = new TLegend(0.1, 0.6, 0.5, 0.9);
    //         gPad->SetLogy();
    //         DataD0JetZPerJet[cent][pT-1] = (TH1D*)AvgJetZPerJet[pT-1][cent]->Clone(Form("DataD0JetZ_%i_%i", cent, pT));
    //         SetName(DataD0JetZPerJet[cent][pT-1], Form("DataD0JetZPerJet_%i_%i", cent, pT));
    //         // DataD0JetZPerJet[cent][pT-1]->Scale(2*TMath::Pi()/taa[cent]);
    //         DataD0JetZPerJet[cent][pT-1]->GetYaxis()->SetRangeUser(1e-4, 1e4);
    //         DataD0JetZPerJet[cent][pT-1]->Draw("E1");
    //         SetColor(TheoryD0JetZPerJet[cent][pT-1], kRed, 20, 1, kRed);
    //         TheoryD0JetZPerJet[cent][pT-1]->Draw("E5 SAME");
    //         theoryvdatazperjetlegend[pT-1][cent]->SetHeader(Form("%i < D0 pT [GeV/#it{c}] < %i", pT, 10), "C");
    //         theoryvdatazperjetlegend[pT-1][cent]->AddEntry(DataD0JetZPerJet[cent][pT-1], "STAR AuAu 0-10%", "P");
    //         theoryvdatazperjetlegend[pT-1][cent]->AddEntry(TheoryD0JetZPerJet[cent][pT-1], "LIDO AuAu 0-10% (MPI = off)", "L");
    //         theoryvdatazperjetlegend[pT-1][cent]->Draw("SAME");
    //         gPad->RedrawAxis();
    //     }
    // }

    // TheoryD0JetZPerJetCanvas->SaveAs(Form("%s/TheoryD0JetZPerJetCanvas.pdf", PlotDir.Data()));

    TLegend *theoryvdatadrperjetlegend[5][3];

    TH1D *Max_Total_dR_Scaled[5][3];
    TH1D *Min_Total_dR_Scaled[5][3];

    TCanvas *TheoryD0JetdRPerJetCanvas = new TCanvas("TheoryD0JetdRPerJetCanvas", "TheoryD0JetdRPerJetCanvas", 2100, 700);
    TheoryD0JetdRPerJetCanvas->Divide(5, 3);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            TheoryD0JetdRPerJetCanvas->cd(cent*5+pT);
            theoryvdatadrperjetlegend[pT-1][cent] = new TLegend(0.1, 0.1, 0.5, 0.4);
            gPad->SetLogy();
            DataD0JetdR[cent][pT-1] = (TH1D*)AvgJetdR[pT-1][cent]->Clone(Form("DataD0JetZ_%i_%i", cent, pT));
            SetName(DataD0JetdR[cent][pT-1], Form("DataD0JetdR_%i_%i", cent, pT));
            // DataD0JetZPerJet[cent][pT-1]->Scale(2*TMath::Pi()/taa[cent]);
            DataD0JetdR[cent][pT-1]->GetXaxis()->SetRangeUser(0, 0.2);
            Max_Total_dR_Scaled[pT-1][cent] = (TH1D*)Max_Total_dR[pT-1][cent]->Clone(Form("Max_Total_dR_Scaled_%i_%i", cent, pT));
            SetName(Max_Total_dR_Scaled[pT-1][cent], Form("Max_Total_dR_Scaled_%i_%i", cent, pT));
            Min_Total_dR_Scaled[pT-1][cent] = (TH1D*)Min_Total_dR[pT-1][cent]->Clone(Form("Min_Total_dR_Scaled_%i_%i", cent, pT));
            SetName(Min_Total_dR_Scaled[pT-1][cent], Form("Min_Total_dR_Scaled_%i_%i", cent, pT));
            Max_Total_dR_Scaled[pT-1][cent]->GetYaxis()->SetRangeUser(1e-4, 1e4);
            Max_Total_dR_Scaled[pT-1][cent]->Draw("HIST ][");
            Min_Total_dR_Scaled[pT-1][cent]->Draw("SAME ][");
            DataD0JetdR[cent][pT-1]->Draw("E1 SAME");
            SetColor(TheoryD0JetdR[cent][pT-1], kRed, 20, 1, kRed);
            TheoryD0JetdR[cent][pT-1]->Draw("E5 SAME");
            theoryvdatadrperjetlegend[pT-1][cent]->SetHeader(Form("%i < D0 pT [GeV/#it{c}] < %i", pT, 10), "C");
            theoryvdatadrperjetlegend[pT-1][cent]->AddEntry(DataD0JetdR[cent][pT-1], "STAR AuAu 0-10%", "P");
            theoryvdatadrperjetlegend[pT-1][cent]->AddEntry(TheoryD0JetdR[cent][pT-1], "LIDO AuAu 0-10% (MPI = off)", "LP");
            theoryvdatadrperjetlegend[pT-1][cent]->Draw("SAME");
            gPad->RedrawAxis();
        }
    }

    TheoryD0JetdRPerJetCanvas->SaveAs(Form("%s/TheoryD0JetdRPerJetCanvas.pdf", PlotDir.Data()));

    TLegend *theoryvdatarcpjetptlegend[5][2];

    TCanvas *TheoryRCPJetPtCanvas = new TCanvas("TheoryRCPJetPtCanvas", "TheoryRCPJetPtCanvas", 2100, 700);
    TheoryRCPJetPtCanvas->Divide(5, 2);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){
            theoryvdatarcpjetptlegend[pT-1][cent] = new TLegend(0.1, 0.6, 0.5, 0.9);
            TheoryRCPJetPtCanvas->cd(cent*5+pT);
            RCP_Max_Total_Pt[pT-1][cent]->Draw("HIST ][");
            RCP_Min_Total_Pt[pT-1][cent]->Draw("SAME ][");
            AvgRCP_JetpT[pT-1][cent]->GetYaxis()->SetRangeUser(0, 2);
            // gPad->SetLogy();
            SetName(TheoryRCP_JetPt[cent][pT-1], Form("RCP_JetPt_%i", pT));
            AvgRCP_JetpT[pT-1][cent]->Draw("E1 SAME");
            TheoryRCP_JetPt[cent][pT-1]->GetYaxis()->SetRangeUser(0,2);
            SetColor(TheoryRCP_JetPt[cent][pT-1], kRed, 20, 1, kRed);
            TheoryRCP_JetPt[cent][pT-1]->Draw("E5 SAME");
            theoryvdatarcpjetptlegend[pT-1][cent]->SetHeader(Form("%i < D0 pT [GeV/#it{c}] < %i", pT, 10), "C");
            theoryvdatarcpjetptlegend[pT-1][cent]->AddEntry(AvgRCP_JetpT[pT-1][cent], "STAR AuAu 0-10%", "P");
            theoryvdatarcpjetptlegend[pT-1][cent]->AddEntry(TheoryRCP_JetPt[cent][pT-1], "LIDO AuAu 0-10% (MPI = off)", "LP");
            theoryvdatarcpjetptlegend[pT-1][cent]->Draw("SAME");
            gPad->RedrawAxis();   
        }
    }

    TheoryRCPJetPtCanvas->SaveAs(Form("%s/TheoryRCPJetPtCanvas.pdf", PlotDir.Data()));


    TH1D *hTheoryD0JetPt[3][5];

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            hTheoryD0JetPt[cent][pT-1] = new TH1D(Form("JetPt_%i_D0pT_%i", cent, pT), Form("JetPt_%i_D0pT_%i", cent, pT), nbins_jpt, binning_jpt);
        }
    }        

    ReadTheoryFiles("LIDO-projection/0-10/D0-jet-spectra-pTD-1-10.dat", hTheoryD0JetPt[0][0], 1.);
    ReadTheoryFiles("LIDO-projection/0-10/D0-jet-spectra-pTD-4-10.dat", hTheoryD0JetPt[0][3], 1.);
    ReadTheoryFiles("LIDO-projection-AuAu200-0-10/D0-jet-spectra-pTD-5-10.dat", hTheoryD0JetPt[0][4], 1.);

    ReadTheoryFiles("LIDO-projection/40-80/D0-jet-spectra-pTD-1-10.dat", hTheoryD0JetPt[2][0], 1.);
    ReadTheoryFiles("LIDO-projection/40-80/D0-jet-spectra-pTD-4-10.dat", hTheoryD0JetPt[2][3], 1.);

    TCanvas *Central = new TCanvas("Central", "Central", 2100, 700);
    Central->cd();
    gPad->SetLogy();

    SetColor(hTheoryD0JetPt[0][0], kRed, 20, 1, kRed);
    hTheoryD0JetPt[0][0]->Draw("E1");
    SetColor(hTheoryD0JetPt[0][3], kBlue, 20, 1, kBlue);
    hTheoryD0JetPt[0][3]->Draw("E1 SAME");
    SetColor(hTheoryD0JetPt[0][4], kGreen, 20, 1, kGreen);
    hTheoryD0JetPt[0][4]->Draw("E1 SAME");
    gPad->BuildLegend();

    TCanvas *Peripheral = new TCanvas("Peripheral", "Peripheral", 2100, 700);
    Peripheral->cd();
    gPad->SetLogy();

    SetColor(hTheoryD0JetPt[2][0], kRed, 20, 1, kRed);
    hTheoryD0JetPt[2][0]->Draw("E1");
    SetColor(hTheoryD0JetPt[2][3], kBlue, 20, 1, kBlue);
    hTheoryD0JetPt[2][3]->Draw("E1 SAME");
    gPad->BuildLegend();

}
