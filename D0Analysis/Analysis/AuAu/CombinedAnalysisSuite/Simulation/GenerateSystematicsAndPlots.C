using namespace std;

#include "BinDef.h"
#include "NewBinDef.h"

#include "iostream"
#include "fstream"

void MakeTmpGraphs(TGraphAsymmErrors *g){
    g = new TGraphAsymmErrors();
    g->SetPoint(0, 1, 0);
}

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

void ReadPPTheoryFiles(TString FileName, TH1D *h, double scalefactor = 1.0){
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
            h->SetBinContent(i-2, pp);
            h->SetBinError(i-2, pperr);
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

// Turns Histogram into a TGraphAsymmErrors
pair <TGraphAsymmErrors *, TGraphAsymmErrors *>GraphPlots(TH1 *Main, TH1 *Center, vector<TH1 *>low, vector<TH1 *>high){
    assert(low.size() == high.size());
    assert(Main->GetNbinsX() == Center->GetNbinsX());

    // cout << "Main: " << Main->GetName() << "\t" << Main->GetNbinsX() << endl;
    // cout << "Center: " << Center->GetName() << "\t" << Center->GetNbinsX() << endl;
    // cout << "Getting Errors From " << low.size() << " sources." << endl;

    TH1D *TotalLow = (TH1D *)low[0]->Clone();
    TotalLow->Reset();
    for (int j = 1; j <= low[0]->GetNbinsX(); j++){
        double totalerr = 0;
        for (int i = 0; i < low.size(); i++){     
            totalerr += pow(low[i]->GetBinContent(j), 2);
        }
        totalerr = sqrt(totalerr);
        // cout << totalerr << endl;
        TotalLow->SetBinContent(j, totalerr);
    }

    TH1D *TotalHigh = (TH1D *)high[0]->Clone();
    TotalHigh->Reset();
    for (int j = 1; j <= high[0]->GetNbinsX(); j++){
        double totalerr = 0;
        for (int i = 0; i < high.size(); i++){     
            totalerr += pow(high[i]->GetBinContent(j), 2);
        }
        totalerr = sqrt(totalerr);
        // cout << totalerr << endl;
        TotalHigh->SetBinContent(j, totalerr);
    }

    TGraphAsymmErrors *stat = new TGraphAsymmErrors(Main);
    stat->SetNameTitle(Form("%s_stat", Main->GetName()), Form("%s_stat", Main->GetName()));
    // cout << stat->GetName() << "\t" << stat->GetN() << endl;

    for (int j = 1; j <= Main->GetNbinsX(); j++){
        // cout << j << "\t" << stat->GetPointX(j-1) << "\t" << stat->GetErrorY(j-1) << endl;
        stat->SetPointEYlow(j-1, Main->GetBinError(j));
        stat->SetPointEYhigh(j-1, Main->GetBinError(j));
    }
    
    TGraphAsymmErrors *sys = new TGraphAsymmErrors(Main);
    sys->SetNameTitle(Form("%s_asymmsys", Main->GetName()), Form("%s_asymmsys", Main->GetName()));
    // cout << sys->GetName() << "\t" << sys->GetN() << endl;

    for (int j = 1; j <= Main->GetNbinsX(); j++){
        // cout << j << "\t" << 0.01*TotalLow->GetBinContent(j)*Main->GetBinContent(j) << endl;
        sys->SetPointEYlow(j-1, 0.01*TotalLow->GetBinContent(j)*Main->GetBinContent(j));
        sys->SetPointEYhigh(j-1, 0.01*TotalHigh->GetBinContent(j)*Main->GetBinContent(j));
    }

    TGraphAsymmErrors *symmstat = new TGraphAsymmErrors(Main);
    symmstat->SetNameTitle(Form("%s_symmstat", Main->GetName()), Form("%s_symmstat", Main->GetName()));
    // cout << stat->GetName() << "\t" << stat->GetN() << endl;

    for (int j = 1; j <= Main->GetNbinsX(); j++){
        // cout << j << "\t" << stat->GetPointX(j-1) << "\t" << stat->GetErrorY(j-1) << endl;
        symmstat->SetPointEYlow(j-1, Main->GetBinError(j));
        symmstat->SetPointEYhigh(j-1, Main->GetBinError(j));
    }

    TGraphAsymmErrors *symmsys = new TGraphAsymmErrors(Main);
    symmsys->SetNameTitle(Form("%s_symmsys", Main->GetName()), Form("%s_symmsys", Main->GetName()));

    // cout << symmsys->GetName() << "\t" << symmsys->GetN() << endl;

    for (int j = 1; j <= Main->GetNbinsX(); j++){
        // cout << j << "\t" << 0.01*TotalLow->GetBinContent(j)*Main->GetBinContent(j) << endl;
        double xlow = Main->GetBinContent(j) - sys->GetErrorYlow(j-1);
        double xhigh = Main->GetBinContent(j) + sys->GetErrorYhigh(j-1);
        cout << xlow << "\t" << xhigh << "\t" << Main->GetBinContent(j) << endl;
        symmsys->SetPoint(j-1, Main->GetBinCenter(j), (xlow+xhigh)/2.);
        symmstat->SetPoint(j-1, Main->GetBinCenter(j), (xlow+xhigh)/2.);
        symmsys->SetPointEYlow(j-1, sys->GetErrorY(j-1));
        symmsys->SetPointEYhigh(j-1, sys->GetErrorY(j-1));
    }

    // cout << symmsys->GetName() << "\t" << symmsys->GetN() << endl;

    pair <TGraphAsymmErrors *, TGraphAsymmErrors *> g;
    g.first = symmstat;
    g.second = symmsys;
    return g;
}


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

int color[12] = {kRed, kGreen-2, kBlue, kBlack, kCyan, kMagenta, kOrange, kOrange+2, kGray+2, kGray+2, kCyan+2, kCyan+2};
int markerstyle[12] = { 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 32};

double taa[3] = {941.23714, 391.35550, 56.62475};
double nevents[3] = {1.0318440e+08, 3.2123506e+08, 4.6679240e+08};

TString JetPtXaxisName = "p_{T,Jet} [GeV/#it{c}]";
TString JetZXaxisName = "z_{Jet} = #vec{p}_{T,Jet}.#vec{p}_{T,D^{0}}/p_{T,Jet}^{2}";
TString JetdRXaxisName = "#Deltar";

TString JetPtYaxisName = "#frac{1}{N_{Evt}} #frac{d^{2}N}{p_{T}dp_{T,Jet}d#eta_{Jet}} [GeV/#it{c}]^{-2}";
TString JetZYaxisName = "#frac{1}{N_{Evt}} #frac{d^{2}N}{zdz_{Jet}d#eta_{Jet}}";
TString JetdRYaxisName = "#frac{1}{N_{Jet}} #frac{dN_{Jet}}{d(#Deltar)}";
TString SystematicUncYaxisName = "(Variation - Nominal)/Nominal";

TString RCPYaxisName = "R_{CP}";
TString RMPYaxisName = "R_{MP}";

TString CentName[4] = {"0-10%", "10-40%", "40-80%", "pp"};

void GenerateSystematicsAndPlots(){
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
    // D0 Reconstruction Efficiency No Vertex Correction Systematic
    AddASet("Aug14_FONLL_D0YieldNoVCEffLower", "Data", 12, 7, 4);
    AddASet("Aug14_FONLL_D0YieldNoVCEffUpper", "Data", 12, 7, 4);
    // D0 Reconstruction Efficiency Vertex Correction Systematic
    AddASet("Aug14_FONLL_D0YieldVCEffLower", "Data", 12, 7, 4);
    AddASet("Aug14_FONLL_D0YieldVCEffUpper", "Data", 12, 7, 4);


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
    
    vector<TString> SystematicName = {"Iter=3", "Iter=5", "Data-Weighed", "PYTHIA-Fit", "Tracking-Efficiency", "D0-Yield-Lower", "D0-Yield-Upper", "D0-Efficiency-No-VC-Lower", "D0-Efficiency-No-VC-Upper", "D0-Efficiency-VC-Lower", "D0-Efficiency-VC-Upper"};
    vector<int> SystematicHist = {3, 4, 1, 2, 5, 6, 7, 8, 9, 10, 11};

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

    // for (int pT = 1; pT <= 5; pT++){
    //     for (int cent = 0; cent < 3; cent++){
    //         AvgRCP_JetpT[pT-1][cent] = (TH1D *)AvgJetpT[pT-1][cent]->Clone();
    //         SetName(AvgRCP_JetpT[pT-1][cent], Form("%s Avg D0pT = %i GeV", RCP_Pt[pT-1][0][cent]->GetName(), pT));
    //         AvgRCP_JetpT[pT-1][cent]->Divide(AvgJetpT[pT-1][2]);
    //         AvgRCP_JetpT[pT-1][cent]->Scale(taa[2]/taa[cent]);

    //         SetBinError(AvgRCP_JetpT[pT-1][cent], RCP_Pt[pT-1][0][cent]);

    //         AvgRCP_JetZ[pT-1][cent] = (TH1D *)AvgJetZ[pT-1][cent]->Clone();
    //         SetName(AvgRCP_JetZ[pT-1][cent], Form("%s Avg D0pT = %i GeV", RCP_Z[pT-1][0][cent]->GetName(), pT));
    //         AvgRCP_JetZ[pT-1][cent]->Divide(AvgJetZ[pT-1][2]);
    //         AvgRCP_JetZ[pT-1][cent]->Scale(taa[2]/taa[cent]);

    //         SetBinError(AvgRCP_JetZ[pT-1][cent], RCP_Z[pT-1][0][cent]);

    //         AvgRCP_JetZPerJet[pT-1][cent] = (TH1D *)AvgJetZPerJet[pT-1][cent]->Clone();
    //         SetName(AvgRCP_JetZPerJet[pT-1][cent], Form("%s Avg D0pT = %i GeV", RCP_ZPerJet[pT-1][0][cent]->GetName(), pT));
    //         AvgRCP_JetZPerJet[pT-1][cent]->Divide(AvgJetZPerJet[pT-1][2]);
    //         AvgRCP_JetZPerJet[pT-1][cent]->Scale(taa[2]/taa[cent]);

    //         SetBinError(AvgRCP_JetZPerJet[pT-1][cent], RCP_ZPerJet[pT-1][0][cent]);

    //         AvgRCP_JetdR[pT-1][cent] = (TH1D *)AvgJetdR[pT-1][cent]->Clone();
    //         SetName(AvgRCP_JetdR[pT-1][cent], Form("%s Avg D0pT = %i GeV", RCP_dR[pT-1][0][cent]->GetName(), pT));
    //         AvgRCP_JetdR[pT-1][cent]->Divide(AvgJetdR[pT-1][2]);

    //         SetBinError(AvgRCP_JetdR[pT-1][cent], RCP_dR[pT-1][0][cent]);
    //     }
    // } 


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

    TCanvas *RCP_PtCanvas_AvgComparison = new TCanvas("RCP_PtCanvas_AvgComparison", "RCP_PtCanvas_AvgComparison", 2100, 700);
    RCP_PtCanvas_AvgComparison->Divide(5,2);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){
            RCP_PtCanvas_AvgComparison->cd(cent*5+pT);
            RCP_Pt[pT-1][0][cent]->Draw("SAME");
            AvgRCP_JetpT[pT-1][cent]->Draw("SAME");
        }
        
    }

    RCP_PtCanvas_AvgComparison->SaveAs(Form("%s/RCP_PtCanvas_AvgComparison_Pt_Iter_%i.pdf", PlotDir.Data(),4));

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

    TCanvas *RCP_ZCanvas_AvgComparison = new TCanvas("RCP_ZCanvas_AvgComparison", "RCP_ZCanvas_AvgComparison", 2100, 700);
    RCP_ZCanvas_AvgComparison->Divide(5,2);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){
            RCP_ZCanvas_AvgComparison->cd(cent*5+pT);
            RCP_Z[pT-1][0][cent]->Draw("SAME");
            AvgRCP_JetZ[pT-1][cent]->Draw("SAME");
        }
        
    }

    RCP_ZCanvas_AvgComparison->SaveAs(Form("%s/RCP_ZCanvas_AvgComparison_Pt_Iter_%i.pdf", PlotDir.Data(),4));

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

    TCanvas *RCP_dRCanvas_AvgComparison = new TCanvas("RCP_dRCanvas_AvgComparison", "RCP_dRCanvas_AvgComparison", 2100, 700);
    RCP_dRCanvas_AvgComparison->Divide(5,2);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){
            RCP_dRCanvas_AvgComparison->cd(cent*5+pT);
            RCP_dR[pT-1][0][cent]->Draw("SAME");
            AvgRCP_JetdR[pT-1][cent]->Draw("SAME");
        }
        
    }

    RCP_dRCanvas_AvgComparison->SaveAs(Form("%s/RCP_dRCanvas_AvgComparison_Pt_Iter_%i.pdf", PlotDir.Data(),4));

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
    TGraphAsymmErrors *FinalJetPt[5][3];
    TGraphAsymmErrors *FinalJetZ[5][3];
    TGraphAsymmErrors *FinalJetZPerJet[5][3];
    TGraphAsymmErrors *FinalJetdR[5][3];

    TGraphAsymmErrors *FinalRCP_Pt[5][2];
    TGraphAsymmErrors *FinalRCP_Z[5][2];
    TGraphAsymmErrors *FinalRCP_ZPerJet[5][3];
    TGraphAsymmErrors *FinalRCP_dR[5][3];

    TGraphAsymmErrors *FinalJetPt_Sys[5][3];
    TGraphAsymmErrors *FinalJetZ_Sys[5][3];
    TGraphAsymmErrors *FinalJetZPerJet_Sys[5][3];
    TGraphAsymmErrors *FinalJetdR_Sys[5][3];

    TGraphAsymmErrors *FinalRCP_Pt_Sys[5][2];
    TGraphAsymmErrors *FinalRCP_Z_Sys[5][2];
    TGraphAsymmErrors *FinalRCP_ZPerJet_Sys[5][3];
    TGraphAsymmErrors *FinalRCP_dR_Sys[5][3];

    vector<TH1 *>lowJetPt[5][3];
    vector<TH1 *>highJetPt[5][3];

    vector<TH1 *>lowJetZ[5][3];
    vector<TH1 *>highJetZ[5][3];

    vector<TH1 *>lowJetZPerJet[5][3];
    vector<TH1 *>highJetZPerJet[5][3];

    vector<TH1 *>lowJetdR[5][3];
    vector<TH1 *>highJetdR[5][3];

    vector<TH1 *>lowRCP_Pt[5][2];
    vector<TH1 *>highRCP_Pt[5][2];

    vector<TH1 *>lowRCP_Z[5][2];
    vector<TH1 *>highRCP_Z[5][2];

    vector<TH1 *>lowRCP_ZPerJet[5][3];
    vector<TH1 *>highRCP_ZPerJet[5][3];

    vector<TH1 *>lowRCP_dR[5][3];
    vector<TH1 *>highRCP_dR[5][3];

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){

            // Iteration Parameters
            lowJetPt[pT-1][cent].push_back((TH1D *)MinSystematicUncertainty_Pt[pT-1][cent]);
            highJetPt[pT-1][cent].push_back((TH1D *)MaxSystematicUncertainty_Pt[pT-1][cent]);

            // Prior Parameters
            lowJetPt[pT-1][cent].push_back((TH1D *)MinSystematicUncertainty_Prior_Pt[pT-1][cent]);
            highJetPt[pT-1][cent].push_back((TH1D *)MaxSystematicUncertainty_Prior_Pt[pT-1][cent]);

            // Tracking Efficiency
            lowJetPt[pT-1][cent].push_back((TH1D *)SystematicUncertainty_Pt[pT-1][4][cent]);
            highJetPt[pT-1][cent].push_back((TH1D *)SystematicUncertainty_Pt[pT-1][4][cent]);

            // D0 Signal Extraction Efficiency
            lowJetPt[pT-1][cent].push_back((TH1D *)SystematicUncertainty_Pt[pT-1][5][cent]);
            highJetPt[pT-1][cent].push_back((TH1D *)SystematicUncertainty_Pt[pT-1][6][cent]);

            // D0 Signal Reconstruction No VC Efficiency
            lowJetPt[pT-1][cent].push_back((TH1D *)SystematicUncertainty_Pt[pT-1][7][cent]);
            highJetPt[pT-1][cent].push_back((TH1D *)SystematicUncertainty_Pt[pT-1][8][cent]);

            // D0 Signal Reconstruction VC Efficiency
            lowJetPt[pT-1][cent].push_back((TH1D *)SystematicUncertainty_Pt[pT-1][9][cent]);
            highJetPt[pT-1][cent].push_back((TH1D *)SystematicUncertainty_Pt[pT-1][10][cent]);

            tie(FinalJetPt[pT-1][cent], FinalJetPt_Sys[pT-1][cent]) = GraphPlots(JetpT[pT-1][0][cent], AvgJetpT[pT-1][cent], lowJetPt[pT-1][cent], highJetPt[pT-1][cent]);
        }
    }

    // for (int i = 1; i <= AvgJetpT[0][0]->GetNbinsX(); i++){
    //     cout << AvgJetpT[0][0]->GetBinContent(i) << "\t" << FinalJetPt[0][0]->GetPointY(i-1) << "\t" << FinalJetPt_Sys[0][0]->GetPointY(i-1) << "\t" << TmpJetPt[0][0]->GetPointY(i-1) << endl;
    // }
    
    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){

            // Iteration Parameters
            lowJetZ[pT-1][cent].push_back((TH1D *)MinSystematicUncertainty_Z[pT-1][cent]);
            highJetZ[pT-1][cent].push_back((TH1D *)MaxSystematicUncertainty_Z[pT-1][cent]);

            // Prior Parameters
            lowJetZ[pT-1][cent].push_back((TH1D *)MinSystematicUncertainty_Prior_Z[pT-1][cent]);
            highJetZ[pT-1][cent].push_back((TH1D *)MaxSystematicUncertainty_Prior_Z[pT-1][cent]);

            // Tracking Efficiency
            lowJetZ[pT-1][cent].push_back((TH1D *)SystematicUncertainty_Z[pT-1][4][cent]);
            highJetZ[pT-1][cent].push_back((TH1D *)SystematicUncertainty_Z[pT-1][4][cent]);

            // D0 Signal Extraction Efficiency
            lowJetZ[pT-1][cent].push_back((TH1D *)SystematicUncertainty_Z[pT-1][5][cent]);
            highJetZ[pT-1][cent].push_back((TH1D *)SystematicUncertainty_Z[pT-1][6][cent]);

            // D0 Signal Reconstruction No VC Efficiency
            lowJetZ[pT-1][cent].push_back((TH1D *)SystematicUncertainty_Z[pT-1][7][cent]);
            highJetZ[pT-1][cent].push_back((TH1D *)SystematicUncertainty_Z[pT-1][8][cent]);

            // D0 Signal Reconstruction VC Efficiency
            lowJetZ[pT-1][cent].push_back((TH1D *)SystematicUncertainty_Z[pT-1][9][cent]);
            highJetZ[pT-1][cent].push_back((TH1D *)SystematicUncertainty_Z[pT-1][10][cent]);

            tie(FinalJetZ[pT-1][cent], FinalJetZ_Sys[pT-1][cent]) = GraphPlots(JetZ[pT-1][0][cent], AvgJetZ[pT-1][cent], lowJetZ[pT-1][cent], highJetZ[pT-1][cent]);
        }
    }

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){

            // Iteration Parameters
            lowJetZPerJet[pT-1][cent].push_back((TH1D *)MinSystematicUncertainty_ZPerJet[pT-1][cent]);
            highJetZPerJet[pT-1][cent].push_back((TH1D *)MaxSystematicUncertainty_ZPerJet[pT-1][cent]);

            // Prior Parameters
            lowJetZPerJet[pT-1][cent].push_back((TH1D *)MinSystematicUncertainty_Prior_ZPerJet[pT-1][cent]);
            highJetZPerJet[pT-1][cent].push_back((TH1D *)MaxSystematicUncertainty_Prior_ZPerJet[pT-1][cent]);

            // Tracking Efficiency
            lowJetZPerJet[pT-1][cent].push_back((TH1D *)SystematicUncertainty_ZPerJet[pT-1][4][cent]);
            highJetZPerJet[pT-1][cent].push_back((TH1D *)SystematicUncertainty_ZPerJet[pT-1][4][cent]);

            // D0 Signal Extraction Efficiency
            lowJetZPerJet[pT-1][cent].push_back((TH1D *)SystematicUncertainty_ZPerJet[pT-1][5][cent]);
            highJetZPerJet[pT-1][cent].push_back((TH1D *)SystematicUncertainty_ZPerJet[pT-1][6][cent]);

            // D0 Signal Reconstruction No VC Efficiency
            lowJetZPerJet[pT-1][cent].push_back((TH1D *)SystematicUncertainty_ZPerJet[pT-1][7][cent]);
            highJetZPerJet[pT-1][cent].push_back((TH1D *)SystematicUncertainty_ZPerJet[pT-1][8][cent]);

            // D0 Signal Reconstruction VC Efficiency
            lowJetZPerJet[pT-1][cent].push_back((TH1D *)SystematicUncertainty_ZPerJet[pT-1][9][cent]);
            highJetZPerJet[pT-1][cent].push_back((TH1D *)SystematicUncertainty_ZPerJet[pT-1][10][cent]);

            tie(FinalJetZPerJet[pT-1][cent], FinalJetZPerJet_Sys[pT-1][cent]) = GraphPlots(JetZPerJet[pT-1][0][cent], AvgJetZPerJet[pT-1][cent], lowJetZPerJet[pT-1][cent], highJetZPerJet[pT-1][cent]);

        }
    }

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){

            // Iteration Parameters
            lowJetdR[pT-1][cent].push_back((TH1D *)MinSystematicUncertainty_dR[pT-1][cent]);
            highJetdR[pT-1][cent].push_back((TH1D *)MaxSystematicUncertainty_dR[pT-1][cent]);

            // Prior Parameters
            lowJetdR[pT-1][cent].push_back((TH1D *)MinSystematicUncertainty_Prior_dR[pT-1][cent]);
            highJetdR[pT-1][cent].push_back((TH1D *)MaxSystematicUncertainty_Prior_dR[pT-1][cent]);

            // Tracking Efficiency
            lowJetdR[pT-1][cent].push_back((TH1D *)SystematicUncertainty_dR[pT-1][4][cent]);
            highJetdR[pT-1][cent].push_back((TH1D *)SystematicUncertainty_dR[pT-1][4][cent]);

            // D0 Signal Extraction Efficiency
            lowJetdR[pT-1][cent].push_back((TH1D *)SystematicUncertainty_dR[pT-1][5][cent]);
            highJetdR[pT-1][cent].push_back((TH1D *)SystematicUncertainty_dR[pT-1][6][cent]);

            // D0 Signal Reconstruction No VC Efficiency
            lowJetdR[pT-1][cent].push_back((TH1D *)SystematicUncertainty_dR[pT-1][7][cent]);
            highJetdR[pT-1][cent].push_back((TH1D *)SystematicUncertainty_dR[pT-1][8][cent]);

            // D0 Signal Reconstruction VC Efficiency
            lowJetdR[pT-1][cent].push_back((TH1D *)SystematicUncertainty_dR[pT-1][9][cent]);
            highJetdR[pT-1][cent].push_back((TH1D *)SystematicUncertainty_dR[pT-1][10][cent]);

            tie(FinalJetdR[pT-1][cent], FinalJetdR_Sys[pT-1][cent]) = GraphPlots(JetdR[pT-1][0][cent], AvgJetdR[pT-1][cent], lowJetdR[pT-1][cent], highJetdR[pT-1][cent]);

        }
    }

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){

            // Iteration Parameters
            lowRCP_Pt[pT-1][cent].push_back((TH1D *)MinSystematicUncertainty_RCP_Pt[pT-1][cent]);
            highRCP_Pt[pT-1][cent].push_back((TH1D *)MaxSystematicUncertainty_RCP_Pt[pT-1][cent]);

            // Prior Parameters
            lowRCP_Pt[pT-1][cent].push_back((TH1D *)MinSystematicUncertainty_Prior_RCP_Pt[pT-1][cent]);
            highRCP_Pt[pT-1][cent].push_back((TH1D *)MaxSystematicUncertainty_Prior_RCP_Pt[pT-1][cent]);

            // Tracking Efficiency
            lowRCP_Pt[pT-1][cent].push_back((TH1D *)SystematicUncertainty_RCP_Pt[pT-1][4][cent]);
            highRCP_Pt[pT-1][cent].push_back((TH1D *)SystematicUncertainty_RCP_Pt[pT-1][4][cent]);

            // D0 Signal Extraction Efficiency (This is uncorrelated, ergo get the error from the jet pT central and peripheral. They will just be added in quadrature.)
            // Central/Midcentral
            lowRCP_Pt[pT-1][cent].push_back((TH1D *)SystematicUncertainty_Pt[pT-1][5][cent]);
            highRCP_Pt[pT-1][cent].push_back((TH1D *)SystematicUncertainty_Pt[pT-1][6][cent]);
            // Peripheral
            lowRCP_Pt[pT-1][cent].push_back((TH1D *)SystematicUncertainty_Pt[pT-1][5][2]);
            highRCP_Pt[pT-1][cent].push_back((TH1D *)SystematicUncertainty_Pt[pT-1][6][2]);

            // D0 Signal Reconstruction No VC Efficiency
            lowRCP_Pt[pT-1][cent].push_back((TH1D *)SystematicUncertainty_RCP_Pt[pT-1][7][cent]);
            highRCP_Pt[pT-1][cent].push_back((TH1D *)SystematicUncertainty_RCP_Pt[pT-1][8][cent]);

            // D0 Signal Reconstruction VC Efficiency (This is uncorrelated, ergo get the error from the jet pT central and peripheral. They will just be added in quadrature.)
            // Central/Midcentral
            lowRCP_Pt[pT-1][cent].push_back((TH1D *)SystematicUncertainty_Pt[pT-1][9][cent]);
            highRCP_Pt[pT-1][cent].push_back((TH1D *)SystematicUncertainty_RCP_Pt[pT-1][10][cent]);
            // Peripheral
            lowRCP_Pt[pT-1][cent].push_back((TH1D *)SystematicUncertainty_Pt[pT-1][9][2]);
            highRCP_Pt[pT-1][cent].push_back((TH1D *)SystematicUncertainty_RCP_Pt[pT-1][10][2]);

            tie(FinalRCP_Pt[pT-1][cent], FinalRCP_Pt_Sys[pT-1][cent]) = GraphPlots(RCP_Pt[pT-1][0][cent], AvgRCP_JetpT[pT-1][cent], lowRCP_Pt[pT-1][cent], highRCP_Pt[pT-1][cent]);
        }
    }

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){

            // Iteration Parameters
            lowRCP_Z[pT-1][cent].push_back((TH1D *)MinSystematicUncertainty_RCP_Z[pT-1][cent]);
            highRCP_Z[pT-1][cent].push_back((TH1D *)MaxSystematicUncertainty_RCP_Z[pT-1][cent]);

            // Prior Parameters
            lowRCP_Z[pT-1][cent].push_back((TH1D *)MinSystematicUncertainty_Prior_RCP_Z[pT-1][cent]);
            highRCP_Z[pT-1][cent].push_back((TH1D *)MaxSystematicUncertainty_Prior_RCP_Z[pT-1][cent]);

            // Tracking Efficiency
            lowRCP_Z[pT-1][cent].push_back((TH1D *)SystematicUncertainty_RCP_Z[pT-1][4][cent]);
            highRCP_Z[pT-1][cent].push_back((TH1D *)SystematicUncertainty_RCP_Z[pT-1][4][cent]);

            // D0 Signal Extraction Efficiency (This is uncorrelated, ergo get the error from the jet pT central and peripheral. They will just be added in quadrature.)
            // Central/Midcentral
            lowRCP_Z[pT-1][cent].push_back((TH1D *)SystematicUncertainty_Z[pT-1][5][cent]);
            highRCP_Pt[pT-1][cent].push_back((TH1D *)SystematicUncertainty_Z[pT-1][6][cent]);
            // Peripheral
            lowRCP_Z[pT-1][cent].push_back((TH1D *)SystematicUncertainty_Z[pT-1][5][2]);
            highRCP_Z[pT-1][cent].push_back((TH1D *)SystematicUncertainty_Z[pT-1][6][2]);

            // D0 Signal Reconstruction No VC Efficiency
            lowRCP_Z[pT-1][cent].push_back((TH1D *)SystematicUncertainty_RCP_Z[pT-1][7][cent]);
            highRCP_Z[pT-1][cent].push_back((TH1D *)SystematicUncertainty_RCP_Z[pT-1][8][cent]);

            // D0 Signal Reconstruction VC Efficiency (This is uncorrelated, ergo get the error from the jet pT central and peripheral. They will just be added in quadrature.)
            // Central/Midcentral
            lowRCP_Z[pT-1][cent].push_back((TH1D *)SystematicUncertainty_Z[pT-1][9][cent]);
            highRCP_Z[pT-1][cent].push_back((TH1D *)SystematicUncertainty_RCP_Z[pT-1][10][cent]);
            // Peripheral
            lowRCP_Z[pT-1][cent].push_back((TH1D *)SystematicUncertainty_Z[pT-1][9][2]);
            highRCP_Z[pT-1][cent].push_back((TH1D *)SystematicUncertainty_RCP_Z[pT-1][10][2]);

            tie(FinalRCP_Z[pT-1][cent], FinalRCP_Z_Sys[pT-1][cent]) = GraphPlots(RCP_Z[pT-1][0][cent], AvgRCP_JetZ[pT-1][cent], lowRCP_Z[pT-1][cent], highRCP_Z[pT-1][cent]);

        }
    }

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){

            // Iteration Parameters
            lowRCP_ZPerJet[pT-1][cent].push_back((TH1D *)MinSystematicUncertainty_RCP_ZPerJet[pT-1][cent]);
            highRCP_ZPerJet[pT-1][cent].push_back((TH1D *)MaxSystematicUncertainty_RCP_ZPerJet[pT-1][cent]);

            // Prior Parameters
            lowRCP_ZPerJet[pT-1][cent].push_back((TH1D *)MinSystematicUncertainty_Prior_RCP_ZPerJet[pT-1][cent]);
            highRCP_ZPerJet[pT-1][cent].push_back((TH1D *)MaxSystematicUncertainty_Prior_RCP_ZPerJet[pT-1][cent]);

            // Tracking Efficiency
            lowRCP_ZPerJet[pT-1][cent].push_back((TH1D *)SystematicUncertainty_RCP_ZPerJet[pT-1][4][cent]);
            highRCP_ZPerJet[pT-1][cent].push_back((TH1D *)SystematicUncertainty_RCP_ZPerJet[pT-1][4][cent]);

            // D0 Signal Extraction Efficiency (This is uncorrelated, ergo get the error from the jet pT central and peripheral. They will just be added in quadrature.)
            // Central/Midcentral
            lowRCP_ZPerJet[pT-1][cent].push_back((TH1D *)SystematicUncertainty_ZPerJet[pT-1][5][cent]);
            highRCP_ZPerJet[pT-1][cent].push_back((TH1D *)SystematicUncertainty_ZPerJet[pT-1][6][cent]);
            // Peripheral
            lowRCP_ZPerJet[pT-1][cent].push_back((TH1D *)SystematicUncertainty_ZPerJet[pT-1][5][2]);
            highRCP_ZPerJet[pT-1][cent].push_back((TH1D *)SystematicUncertainty_ZPerJet[pT-1][6][2]);
            
            // D0 Signal Reconstruction No VC Efficiency
            lowRCP_ZPerJet[pT-1][cent].push_back((TH1D *)SystematicUncertainty_RCP_ZPerJet[pT-1][7][cent]);
            highRCP_ZPerJet[pT-1][cent].push_back((TH1D *)SystematicUncertainty_RCP_ZPerJet[pT-1][8][cent]);

            // D0 Signal Reconstruction VC Efficiency (This is uncorrelated, ergo get the error from the jet pT central and peripheral. They will just be added in quadrature.)
            // Central/Midcentral
            lowRCP_ZPerJet[pT-1][cent].push_back((TH1D *)SystematicUncertainty_ZPerJet[pT-1][9][cent]);
            highRCP_ZPerJet[pT-1][cent].push_back((TH1D *)SystematicUncertainty_RCP_ZPerJet[pT-1][10][cent]);
            // Peripheral
            lowRCP_ZPerJet[pT-1][cent].push_back((TH1D *)SystematicUncertainty_ZPerJet[pT-1][9][2]);
            highRCP_ZPerJet[pT-1][cent].push_back((TH1D *)SystematicUncertainty_RCP_ZPerJet[pT-1][10][2]);

            tie(FinalRCP_ZPerJet[pT-1][cent], FinalRCP_ZPerJet_Sys[pT-1][cent]) = GraphPlots(RCP_ZPerJet[pT-1][0][cent], AvgRCP_JetZPerJet[pT-1][cent], lowRCP_ZPerJet[pT-1][cent], highRCP_ZPerJet[pT-1][cent]);

        }
    }

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){

            // Iteration Parameters
            lowRCP_dR[pT-1][cent].push_back((TH1D *)MinSystematicUncertainty_RCP_dR[pT-1][cent]);
            highRCP_dR[pT-1][cent].push_back((TH1D *)MaxSystematicUncertainty_RCP_dR[pT-1][cent]);

            // Prior Parameters
            lowRCP_dR[pT-1][cent].push_back((TH1D *)MinSystematicUncertainty_Prior_RCP_dR[pT-1][cent]);
            highRCP_dR[pT-1][cent].push_back((TH1D *)MaxSystematicUncertainty_Prior_RCP_dR[pT-1][cent]);

            // Tracking Efficiency
            lowRCP_dR[pT-1][cent].push_back((TH1D *)SystematicUncertainty_RCP_dR[pT-1][4][cent]);
            highRCP_dR[pT-1][cent].push_back((TH1D *)SystematicUncertainty_RCP_dR[pT-1][4][cent]);

            // D0 Signal Extraction Efficiency (This is uncorrelated, ergo get the error from the jet pT central and peripheral. They will just be added in quadrature.)
            // Central/Midcentral
            lowRCP_dR[pT-1][cent].push_back((TH1D *)SystematicUncertainty_dR[pT-1][5][cent]);
            highRCP_dR[pT-1][cent].push_back((TH1D *)SystematicUncertainty_dR[pT-1][6][cent]);
            // Peripheral
            lowRCP_dR[pT-1][cent].push_back((TH1D *)SystematicUncertainty_dR[pT-1][5][2]);
            highRCP_dR[pT-1][cent].push_back((TH1D *)SystematicUncertainty_dR[pT-1][6][2]);
            
            // D0 Signal Reconstruction No VC Efficiency
            lowRCP_dR[pT-1][cent].push_back((TH1D *)SystematicUncertainty_RCP_dR[pT-1][7][cent]);
            highRCP_dR[pT-1][cent].push_back((TH1D *)SystematicUncertainty_RCP_dR[pT-1][8][cent]);

            // D0 Signal Reconstruction VC Efficiency (This is uncorrelated, ergo get the error from the jet pT central and peripheral. They will just be added in quadrature.)
            // Central/Midcentral
            lowRCP_dR[pT-1][cent].push_back((TH1D *)SystematicUncertainty_dR[pT-1][9][cent]);
            highRCP_dR[pT-1][cent].push_back((TH1D *)SystematicUncertainty_RCP_dR[pT-1][10][cent]);
            // Peripheral
            lowRCP_dR[pT-1][cent].push_back((TH1D *)SystematicUncertainty_dR[pT-1][9][2]);
            highRCP_dR[pT-1][cent].push_back((TH1D *)SystematicUncertainty_RCP_dR[pT-1][10][2]);

            tie(FinalRCP_dR[pT-1][cent], FinalRCP_dR_Sys[pT-1][cent]) = GraphPlots(RCP_dR[pT-1][0][cent], AvgRCP_JetdR[pT-1][cent], lowRCP_dR[pT-1][cent], highRCP_dR[pT-1][cent]);

        }
    }

    TCanvas *FinalJetPtCanvas = new TCanvas("FinalJetPtCanvas", "FinalJetPtCanvas", 2100, 700);
    FinalJetPtCanvas->Divide(5, 3);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            FinalJetPtCanvas->cd(cent*5+pT);
            gPad->SetLogy();
            // FinalJetPt[pT-1][cent]->SetRange(5, 20);
            FinalJetPt_Sys[pT-1][cent]->SetLineColor(kBlue);
            FinalJetPt_Sys[pT-1][cent]->SetFillColor(kBlue);
            FinalJetPt[pT-1][cent]->SetMarkerStyle(29);
            FinalJetPt[pT-1][cent]->SetMarkerSize(1.2);
            FinalJetPt_Sys[pT-1][cent]->SetFillStyle(0);
            FinalJetPt_Sys[pT-1][cent]->Draw("A2");
            FinalJetPt[pT-1][cent]->Draw("P SAME");
            FinalJetPt_Sys[pT-1][cent]->GetXaxis()->SetRangeUser(5, 20);
            FinalJetPt_Sys[pT-1][cent]->GetYaxis()->SetRangeUser(1e-10,1e-2);
            FinalJetPt_Sys[pT-1][cent]->GetXaxis()->SetTitle(JetPtXaxisName.Data());
            FinalJetPt_Sys[pT-1][cent]->GetYaxis()->SetTitle(JetPtYaxisName.Data());
        }
        
    }
    
    FinalJetPtCanvas->SaveAs(Form("%s/FinalJetPt_Iter_%i.pdf", PlotDir.Data(),4));

    TCanvas *FinalJetZCanvas = new TCanvas("FinalJetZCanvas", "FinalJetZCanvas", 2100, 700);
    FinalJetZCanvas->Divide(5, 3);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            FinalJetZCanvas->cd(cent*5+pT);
            gPad->SetLogy();
            // FinalJetZ[pT-1][cent]->SetRange(5, 20);
            FinalJetZ_Sys[pT-1][cent]->SetLineColor(kBlue);
            FinalJetZ_Sys[pT-1][cent]->SetFillColor(kBlue);
            FinalJetZ[pT-1][cent]->SetMarkerStyle(29);
            FinalJetZ[pT-1][cent]->SetMarkerSize(1.2);
            FinalJetZ_Sys[pT-1][cent]->SetFillStyle(0);
            FinalJetZ_Sys[pT-1][cent]->Draw("A2");
            FinalJetZ[pT-1][cent]->Draw("P SAME");
            FinalJetZ_Sys[pT-1][cent]->GetXaxis()->SetRangeUser(0,1);
            FinalJetZ_Sys[pT-1][cent]->GetYaxis()->SetRangeUser(1e-8,1e2);
            FinalJetZ_Sys[pT-1][cent]->GetXaxis()->SetTitle(JetZXaxisName.Data());
            FinalJetZ_Sys[pT-1][cent]->GetYaxis()->SetTitle(JetZYaxisName.Data());
        }
        
    }

    FinalJetZCanvas->SaveAs(Form("%s/FinalJetZ_Iter_%i.pdf", PlotDir.Data(),4));

    TCanvas *FinalJetdRCanvas = new TCanvas("FinalJetdRCanvas", "FinalJetdRCanvas", 2100, 700);
    FinalJetdRCanvas->Divide(5, 3);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            FinalJetdRCanvas->cd(cent*5+pT);
            gPad->SetLogy();
            // FinalJetdR[pT-1][cent]->SetRange(5, 20);
            FinalJetdR_Sys[pT-1][cent]->SetLineColor(kBlue);
            FinalJetdR_Sys[pT-1][cent]->SetFillColor(kBlue);
            FinalJetdR[pT-1][cent]->SetMarkerStyle(29);
            FinalJetdR[pT-1][cent]->SetMarkerSize(1.2);
            FinalJetdR_Sys[pT-1][cent]->SetFillStyle(0);
            FinalJetdR_Sys[pT-1][cent]->Draw("A2");
            FinalJetdR[pT-1][cent]->Draw("P SAME");
            FinalJetdR_Sys[pT-1][cent]->GetXaxis()->SetRangeUser(0, 0.2);
            FinalJetdR_Sys[pT-1][cent]->GetYaxis()->SetRangeUser(1e-4,1e4);
            FinalJetdR_Sys[pT-1][cent]->GetXaxis()->SetTitle(JetdRXaxisName.Data());
            FinalJetdR_Sys[pT-1][cent]->GetYaxis()->SetTitle(JetdRYaxisName.Data());
        }
        
    }

    FinalJetdRCanvas->SaveAs(Form("%s/FinalJetdR_Iter_%i.pdf", PlotDir.Data(),4));

    TLine *PtLineAtOne = new TLine(5, 1, 20, 1);
    PtLineAtOne->SetLineColor(kBlack);

    TCanvas *FinalRCP_PtCanvas = new TCanvas("FinalRCP_PtCanvas", "FinalRCP_PtCanvas", 2100, 700);
    FinalRCP_PtCanvas->Divide(5, 2);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){
            FinalRCP_PtCanvas->cd(cent*5+pT);
            // gPad->SetLogy();
            // FinalRCP_Pt[pT-1][cent]->SetRange(5, 20);
            FinalRCP_Pt_Sys[pT-1][cent]->SetLineColor(kBlue);
            FinalRCP_Pt_Sys[pT-1][cent]->SetFillColor(kBlue);
            FinalRCP_Pt[pT-1][cent]->SetMarkerStyle(29);
            FinalRCP_Pt[pT-1][cent]->SetMarkerSize(1.2);
            FinalRCP_Pt_Sys[pT-1][cent]->SetFillStyle(0);
            FinalRCP_Pt_Sys[pT-1][cent]->Draw("A2");
            FinalRCP_Pt[pT-1][cent]->Draw("P SAME");
            FinalRCP_Pt_Sys[pT-1][cent]->GetXaxis()->SetRangeUser(5, 20);
            FinalRCP_Pt_Sys[pT-1][cent]->GetYaxis()->SetRangeUser(0.1, 2);
            FinalRCP_Pt_Sys[pT-1][cent]->GetXaxis()->SetTitle(JetPtXaxisName.Data());
            FinalRCP_Pt_Sys[pT-1][cent]->GetYaxis()->SetTitle(RCPYaxisName.Data());
            PtLineAtOne->Draw("SAME");
        }
        
    }

    FinalRCP_PtCanvas->SaveAs(Form("%s/FinalRCP_Pt_Iter_%i.pdf", PlotDir.Data(),4));

    TLine *ZLineAtOne = new TLine(0, 1, 1, 1);
    ZLineAtOne->SetLineColor(kBlack);

    TCanvas *FinalRCP_ZCanvas = new TCanvas("FinalRCP_ZCanvas", "FinalRCP_ZCanvas", 2100, 700);
    FinalRCP_ZCanvas->Divide(5, 2);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){
            FinalRCP_ZCanvas->cd(cent*5+pT);
            // gPad->SetLogy();
            // FinalRCP_Z[pT-1][cent]->SetRange(5, 20);
            FinalRCP_Z_Sys[pT-1][cent]->SetLineColor(kBlue);
            FinalRCP_Z_Sys[pT-1][cent]->SetFillColor(kBlue);
            FinalRCP_Z[pT-1][cent]->SetMarkerStyle(29);
            FinalRCP_Z[pT-1][cent]->SetMarkerSize(1.2);
            FinalRCP_Z_Sys[pT-1][cent]->SetFillStyle(0);
            FinalRCP_Z_Sys[pT-1][cent]->Draw("A2");
            FinalRCP_Z[pT-1][cent]->Draw("P SAME");
            FinalRCP_Z_Sys[pT-1][cent]->GetXaxis()->SetRangeUser(0,1);
            FinalRCP_Z_Sys[pT-1][cent]->GetYaxis()->SetRangeUser(0.1, 2.5);
            FinalRCP_Z_Sys[pT-1][cent]->GetXaxis()->SetTitle(JetZXaxisName.Data());
            FinalRCP_Z_Sys[pT-1][cent]->GetYaxis()->SetTitle(RCPYaxisName.Data());
            ZLineAtOne->Draw("SAME");
        }
        
    }

    FinalRCP_ZCanvas->SaveAs(Form("%s/FinalRCP_Z_Iter_%i.pdf", PlotDir.Data(),4));

    TLine *dRLineAtOne = new TLine(0, 1, 0.2, 1);
    dRLineAtOne->SetLineColor(kBlack);

    TCanvas *FinalRCP_dRCanvas = new TCanvas("FinalRCP_dRCanvas", "FinalRCP_dRCanvas", 2100, 700);
    FinalRCP_dRCanvas->Divide(5, 2);

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){
            FinalRCP_dRCanvas->cd(cent*5+pT);
            // gPad->SetLogy();
            // FinalRCP_dR[pT-1][cent]->SetRange(5, 20);
            FinalRCP_dR_Sys[pT-1][cent]->SetLineColor(kBlue);
            FinalRCP_dR_Sys[pT-1][cent]->SetFillColor(kBlue);
            FinalRCP_dR[pT-1][cent]->SetMarkerStyle(29);
            FinalRCP_dR[pT-1][cent]->SetMarkerSize(1.2);
            FinalRCP_dR_Sys[pT-1][cent]->SetFillStyle(0);
            FinalRCP_dR_Sys[pT-1][cent]->Draw("A2");
            FinalRCP_dR[pT-1][cent]->Draw("P SAME");
            FinalRCP_dR_Sys[pT-1][cent]->GetXaxis()->SetRangeUser(0, 0.2);
            FinalRCP_dR_Sys[pT-1][cent]->GetYaxis()->SetRangeUser(0.1, 2);
            FinalRCP_dR_Sys[pT-1][cent]->GetXaxis()->SetTitle(JetdRXaxisName.Data());
            FinalRCP_dR_Sys[pT-1][cent]->GetYaxis()->SetTitle(RCPYaxisName.Data());
            dRLineAtOne->Draw("SAME");
        }
        
    }

    FinalRCP_dRCanvas->SaveAs(Form("%s/FinalRCP_dR_Iter_%i.pdf", PlotDir.Data(),4));

    TH1D *TheoryD0JetPt[5][4];
    TH1D *TheoryD0JetZ[5][4];
    TH1D *TheoryD0JetZPerJet[5][4];
    TH1D *TheoryD0JetdR[5][4];

    TH1D *TheoryRCP_JetPt[5][2];
    TH1D *TheoryRCP_JetZ[5][2];
    TH1D *TheoryRCP_JetdR[5][2];

    int nz_bins_theory = 9;
    double nbinsz_theory[10] = {0., 0.1, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 4; cent++){
            TheoryD0JetPt[pT-1][cent] = new TH1D(Form("TheoryD0JetPt_Pt_%i_Cent_%i", pT, cent), Form("TheoryD0JetPt_Pt_%i_Cent_%i", pT, cent), nbins_jpt, binning_jpt);
            TheoryD0JetZ[pT-1][cent] = new TH1D(Form("TheoryD0JetZ_Pt_%i_Cent_%i", pT, cent), Form("TheoryD0JetZ_Pt_%i_Cent_%i", pT, cent), nz_bins_theory, nbinsz_theory);
            TheoryD0JetZPerJet[pT-1][cent] = new TH1D(Form("TheoryD0JetZPerJet_Pt_%i_Cent_%i", pT, cent), Form("TheoryD0JetZPerJet_Pt_%i_Cent_%i", pT, cent), nz_bins_theory, nbinsz_theory);
            TheoryD0JetdR[pT-1][cent] = new TH1D(Form("TheoryD0JetdR_Pt_%i_Cent_%i", pT, cent), Form("TheoryD0JetdR_Pt_%i_Cent_%i", pT, cent), ndrbins, drbins);
        }

        for (int cent = 0; cent < 2; cent++){
            TheoryRCP_JetPt[pT-1][cent] = new TH1D(Form("TheoryRCP_JetPt_Pt_%i_Cent_%i", pT, cent), Form("TheoryRCP_JetPt_Pt_%i_Cent_%i", pT, cent), nbins_jpt, binning_jpt);
            TheoryRCP_JetZ[pT-1][cent] = new TH1D(Form("TheoryRCP_JetZ_Pt_%i_Cent_%i", pT, cent), Form("TheoryRCP_JetZ_Pt_%i_Cent_%i", pT, cent), nz_bins_theory, nbinsz_theory);
            TheoryRCP_JetdR[pT-1][cent] = new TH1D(Form("TheoryRCP_JetdR_Pt_%i_Cent_%i", pT, cent), Form("TheoryRCP_JetdR_Pt_%i_Cent_%i", pT, cent), ndrbins, drbins);
        }
    }

    double inelcrosssection = 43.82; // mb

    ReadTheoryFiles("LIDO-projection/0-10/D0-jet-spectra-pTD-1-10.dat", TheoryD0JetPt[0][0], 1./(2.*TMath::Pi()*inelcrosssection));
    ReadTheoryFiles("LIDO-projection/0-10/D0-jet-spectra-pTD-4-10.dat", TheoryD0JetPt[3][0], 1./(2.*TMath::Pi()*inelcrosssection));
    ReadTheoryFiles("LIDO-projection-AuAu200-0-10/D0-jet-spectra-pTD-5-10.dat", TheoryD0JetPt[4][0], 1./(2.*TMath::Pi()*inelcrosssection));

    ReadTheoryFiles("LIDO-projection/40-80/D0-jet-spectra-pTD-1-10.dat", TheoryD0JetPt[0][2], 1./(2.*TMath::Pi()*inelcrosssection));
    ReadTheoryFiles("LIDO-projection/40-80/D0-jet-spectra-pTD-4-10.dat", TheoryD0JetPt[3][2], 1./(2.*TMath::Pi()*inelcrosssection));

    ReadPPTheoryFiles("LIDO-projection/0-10/D0-jet-spectra-pTD-1-10.dat", TheoryD0JetPt[0][3], 1./(2.*TMath::Pi()*inelcrosssection));

    ReadTheoryFiles("LIDO-projection/0-10/D0inJet_z_cross-section-pTD-1-10.dat", TheoryD0JetZ[0][0], 1./(2.*TMath::Pi()*inelcrosssection));
    ReadTheoryFiles("LIDO-projection/0-10/D0inJet_z_cross-section-pTD-4-10.dat", TheoryD0JetZ[3][0], 1./(2.*TMath::Pi()*inelcrosssection));
    ReadTheoryFiles("LIDO-projection-AuAu200-0-10/D0inJet_z_cross-section-pTD-5-10.dat", TheoryD0JetZ[4][0], 1./(2.*TMath::Pi()*inelcrosssection));

    cout << "<==================== D0 Z From LIDO --> 40 - 80 ====================> " << endl;
    ReadTheoryFiles("LIDO-projection/40-80/D0inJet_z_cross-section-pTD-1-10.dat", TheoryD0JetZ[0][2], 1./(2.*TMath::Pi()*inelcrosssection));
    cout << "<==================== D0 Z From LIDO --> 40 - 80 ====================> " << endl;
    ReadTheoryFiles("LIDO-projection/40-80/D0inJet_z_cross-section-pTD-4-10.dat", TheoryD0JetZ[3][2], 1./(2.*TMath::Pi()*inelcrosssection));

    ReadPPTheoryFiles("LIDO-projection/0-10/D0inJet_z_cross-section-pTD-1-10.dat", TheoryD0JetZ[0][3], 1./(2.*TMath::Pi()*inelcrosssection));

    ReadTheoryFiles("LIDO-projection/0-10/D0inJet_z_per-jet-pTD-1-10.dat", TheoryD0JetZPerJet[0][0], 1.);
    ReadTheoryFiles("LIDO-projection/0-10/D0inJet_z_per-jet-pTD-4-10.dat", TheoryD0JetZPerJet[3][0], 1.);

    ReadTheoryFiles("LIDO-projection/40-80/D0inJet_z_per-jet-pTD-1-10.dat", TheoryD0JetZPerJet[0][2], 1.);
    ReadTheoryFiles("LIDO-projection/40-80/D0inJet_z_per-jet-pTD-4-10.dat", TheoryD0JetZPerJet[3][2], 1.);

    cout << "<==================== D0 dR From LIDO --> 0 - 10 ====================> " << endl;
    ReadTheoryFiles("LIDO-projection/0-10/D0inJet_r_per-jet-pTD-1-10.dat", TheoryD0JetdR[0][0], 1.);
    ReadTheoryFiles("LIDO-projection/0-10/D0inJet_r_per-jet-pTD-4-10.dat", TheoryD0JetdR[3][0], 1.);
    // ReadTheoryFiles("LIDO-projection-AuAu200-0-10/D0inJet_r_per-jet-pTD-5-10.dat", TheoryD0JetdR[4][0], 1.);

    cout << "<==================== D0 dR From LIDO --> 40 - 80 ====================> " << endl;
    ReadTheoryFiles("LIDO-projection/40-80/D0inJet_r_per-jet-pTD-1-10.dat", TheoryD0JetdR[0][2], 1.);
    ReadTheoryFiles("LIDO-projection/40-80/D0inJet_r_per-jet-pTD-4-10.dat", TheoryD0JetdR[3][2], 1.);

    ReadPPTheoryFiles("LIDO-projection/0-10/D0inJet_r_per-jet-pTD-1-10.dat", TheoryD0JetdR[0][3], 1.);


    ReadRCPTheoryFiles("LIDO-projection/RCP-1-10.dat", TheoryRCP_JetPt[0][0], 1.);
    ReadRCPTheoryFiles("LIDO-projection/RCP-4-10.dat", TheoryRCP_JetPt[3][0], 1.);

    TheoryRCP_JetZ[0][0] = (TH1D *)TheoryD0JetZ[0][0]->Clone("TheoryRCP_JetZ_Pt_1_Cent_0");
    SetName(TheoryRCP_JetZ[0][0], "TheoryRCP_JetZ_Pt_1_Cent_0");
    TheoryRCP_JetZ[0][0]->Divide(TheoryD0JetZ[0][2]);
    // TheoryRCP_JetZ[0][0]->Scale(taa[2]/taa[0]);
    TheoryRCP_JetZ[3][0] = (TH1D *)TheoryD0JetZ[3][0]->Clone("TheoryRCP_JetZ_Pt_4_Cent_0");
    SetName(TheoryRCP_JetZ[3][0], "TheoryRCP_JetZ_Pt_4_Cent_0");
    TheoryRCP_JetZ[3][0]->Divide(TheoryD0JetZ[3][2]);

    TheoryRCP_JetdR[0][0] = (TH1D *)TheoryD0JetdR[0][0]->Clone("TheoryRCP_JetdR_Pt_1_Cent_0");
    SetName(TheoryRCP_JetdR[0][0], "TheoryRCP_JetdR_Pt_1_Cent_0");
    TheoryRCP_JetdR[0][0]->Divide(TheoryD0JetdR[0][2]);
    TheoryRCP_JetdR[3][0] = (TH1D *)TheoryD0JetdR[3][0]->Clone("TheoryRCP_JetdR_Pt_4_Cent_0");
    SetName(TheoryRCP_JetdR[3][0], "TheoryRCP_JetdR_Pt_4_Cent_0");
    TheoryRCP_JetdR[3][0]->Divide(TheoryD0JetdR[3][2]);

    TCanvas *FinalJetPtCanvas_LIDO = new TCanvas("FinalJetPtCanvas_LIDO", "FinalJetPtCanvas_LIDO", 2100, 700);
    FinalJetPtCanvas_LIDO->Divide(5, 3);

    TLegend *theoryvdatajetptlegend;
    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            FinalJetPtCanvas_LIDO->cd(cent*5+pT);
            if (pT == 1 && cent == 0)theoryvdatajetptlegend = new TLegend(0.4, 0.6, 0.9, 0.9);
            gPad->SetLogy();
            // FinalJetPt[pT-1][cent]->SetRange(5, 20);
            FinalJetPt_Sys[pT-1][cent]->SetLineColor(kBlue);
            FinalJetPt_Sys[pT-1][cent]->SetFillColor(kBlue);
            FinalJetPt[pT-1][cent]->SetMarkerStyle(29);
            FinalJetPt[pT-1][cent]->SetMarkerSize(1.2);
            FinalJetPt_Sys[pT-1][cent]->SetFillStyle(0);
            FinalJetPt_Sys[pT-1][cent]->Draw("A2");
            FinalJetPt[pT-1][cent]->Draw("P SAME");
            FinalJetPt_Sys[pT-1][cent]->GetXaxis()->SetRangeUser(5, 20);
            FinalJetPt_Sys[pT-1][cent]->GetYaxis()->SetRangeUser(1e-10,1e-2);
            TheoryD0JetPt[pT-1][cent]->SetLineColor(kGreen-2);
            TheoryD0JetPt[pT-1][cent]->SetFillColor(kGreen-2);
            TheoryD0JetPt[pT-1][cent]->SetMarkerColor(kGreen-2);
            TheoryD0JetPt[pT-1][cent]->Scale(taa[cent]);
            TheoryD0JetPt[pT-1][cent]->SetMarkerStyle(20);
            TheoryD0JetPt[pT-1][cent]->Draw("E5 SAME");
            if (pT == 1 && cent == 0){
                theoryvdatajetptlegend->AddEntry(FinalJetPt[pT-1][cent], "STAR AuAu", "P");
                theoryvdatajetptlegend->AddEntry(TheoryD0JetPt[pT-1][cent], "LIDO AuAu (MPI = off)", "P");
            }
            theoryvdatajetptlegend->Draw("SAME");
        }
        
    }
    
    FinalJetPtCanvas_LIDO->SaveAs(Form("%s/FinalJetPt_LIDO_Iter_%i.pdf", PlotDir.Data(),4));

    TCanvas *FinalJetZCanvas_LIDO = new TCanvas("FinalJetZCanvas_LIDO", "FinalJetZCanvas_LIDO", 2100, 700);
    FinalJetZCanvas_LIDO->Divide(5, 3);

    TLegend *theoryvdatazlegend;
    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            FinalJetZCanvas_LIDO->cd(cent*5+pT);
            gPad->SetLogy();
            if (pT == 1 && cent == 0) theoryvdatazlegend = new TLegend(0.1, 0.6, 0.5, 0.9);
            // FinalJetZ[pT-1][cent]->SetRange(5, 20);
            FinalJetZ_Sys[pT-1][cent]->SetLineColor(kBlue);
            FinalJetZ_Sys[pT-1][cent]->SetFillColor(kBlue);
            FinalJetZ[pT-1][cent]->SetMarkerStyle(29);
            FinalJetZ[pT-1][cent]->SetMarkerSize(1.2);
            FinalJetZ_Sys[pT-1][cent]->SetFillStyle(0);
            FinalJetZ_Sys[pT-1][cent]->Draw("A2");
            FinalJetZ[pT-1][cent]->Draw("P SAME");
            FinalJetZ_Sys[pT-1][cent]->GetXaxis()->SetRangeUser(0,1);
            FinalJetZ_Sys[pT-1][cent]->GetYaxis()->SetRangeUser(1e-8,1e2);
            TheoryD0JetZ[pT-1][cent]->SetLineColor(kGreen-2);
            TheoryD0JetZ[pT-1][cent]->SetFillColor(kGreen-2);
            TheoryD0JetZ[pT-1][cent]->SetMarkerColor(kGreen-2);
            TheoryD0JetZ[pT-1][cent]->Scale(taa[cent]);
            TheoryD0JetZ[pT-1][cent]->SetMarkerStyle(20);
            TheoryD0JetZ[pT-1][cent]->Draw("E5 SAME");
            if (pT == 1 && cent == 0){
                theoryvdatazlegend->AddEntry(FinalJetZ[pT-1][cent], "STAR AuAu", "P");
                theoryvdatazlegend->AddEntry(TheoryD0JetZ[pT-1][cent], "LIDO AuAu (MPI = off)", "P");
            }
            theoryvdatazlegend->Draw("SAME");
        }
        
    }

    FinalJetZCanvas_LIDO->SaveAs(Form("%s/FinalJetZ_LIDO_Iter_%i.pdf", PlotDir.Data(),4));

    TCanvas *FinalJetdRCanvas_LIDO = new TCanvas("FinalJetdRCanvas_LIDO", "FinalJetdRCanvas_LIDO", 2100, 700);
    FinalJetdRCanvas_LIDO->Divide(5, 3);

    TLegend *theoryvdatadRlegend;
    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            FinalJetdRCanvas_LIDO->cd(cent*5+pT);
            gPad->SetLogy();
            if (pT == 1 && cent == 0) theoryvdatadRlegend = new TLegend(0.1, 0.1, 0.5, 0.4);
            // FinalJetdR[pT-1][cent]->SetRange(5, 20);
            FinalJetdR_Sys[pT-1][cent]->SetLineColor(kBlue);
            FinalJetdR_Sys[pT-1][cent]->SetFillColor(kBlue);
            FinalJetdR[pT-1][cent]->SetMarkerStyle(29);
            FinalJetdR[pT-1][cent]->SetMarkerSize(1.2);
            FinalJetdR_Sys[pT-1][cent]->SetFillStyle(0);
            FinalJetdR_Sys[pT-1][cent]->Draw("A2");
            FinalJetdR[pT-1][cent]->Draw("P SAME");
            FinalJetdR_Sys[pT-1][cent]->GetXaxis()->SetRangeUser(0, 0.2);
            FinalJetdR_Sys[pT-1][cent]->GetYaxis()->SetRangeUser(1e-4,1e4);
            TheoryD0JetdR[pT-1][cent]->SetLineColor(kGreen-2);
            TheoryD0JetdR[pT-1][cent]->SetFillColor(kGreen-2);
            TheoryD0JetdR[pT-1][cent]->SetMarkerColor(kGreen-2);
            TheoryD0JetdR[pT-1][cent]->SetMarkerStyle(20);
            TheoryD0JetdR[pT-1][cent]->Draw("E5 SAME");
            if (pT == 1 && cent == 0){
                theoryvdatadRlegend->AddEntry(FinalJetdR[pT-1][cent], "STAR AuAu", "P");
                theoryvdatadRlegend->AddEntry(TheoryD0JetdR[pT-1][cent], "LIDO AuAu (MPI = off)", "P");
            }
            theoryvdatadRlegend->Draw("SAME");
        }
        
    }

    FinalJetdRCanvas_LIDO->SaveAs(Form("%s/FinalJetdR_LIDO_Iter_%i.pdf", PlotDir.Data(),4));

    TCanvas *FinalRCP_PtCanvas_LIDO = new TCanvas("FinalRCP_PtCanvas_LIDO", "FinalRCP_PtCanvas_LIDO", 2100, 700);
    FinalRCP_PtCanvas_LIDO->Divide(5, 2);

    TLegend *theoryvdatarcpjetptlegend;
    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){
            FinalRCP_PtCanvas_LIDO->cd(cent*5+pT);
            if (pT == 1 && cent == 0) theoryvdatarcpjetptlegend = new TLegend(0.1, 0.6, 0.5, 0.9);
            // gPad->SetLogy();
            // FinalRCP_Pt[pT-1][cent]->SetRange(5, 20);
            FinalRCP_Pt_Sys[pT-1][cent]->SetLineColor(kBlue);
            FinalRCP_Pt_Sys[pT-1][cent]->SetFillColor(kBlue);
            FinalRCP_Pt[pT-1][cent]->SetMarkerStyle(29);
            FinalRCP_Pt[pT-1][cent]->SetMarkerSize(1.2);
            FinalRCP_Pt_Sys[pT-1][cent]->SetFillStyle(0);
            FinalRCP_Pt_Sys[pT-1][cent]->Draw("A2");
            FinalRCP_Pt[pT-1][cent]->Draw("P SAME");
            FinalRCP_Pt_Sys[pT-1][cent]->GetXaxis()->SetRangeUser(5, 20);
            FinalRCP_Pt_Sys[pT-1][cent]->GetYaxis()->SetRangeUser(0.1, 2);
            TheoryRCP_JetPt[pT-1][cent]->SetLineColor(kGreen-2);
            TheoryRCP_JetPt[pT-1][cent]->SetFillColor(kGreen-2);
            TheoryRCP_JetPt[pT-1][cent]->SetMarkerColor(kGreen-2);
            TheoryRCP_JetPt[pT-1][cent]->SetMarkerStyle(20);
            TheoryRCP_JetPt[pT-1][cent]->Draw("E5 SAME");
            PtLineAtOne->Draw("SAME");
            if (pT == 1 && cent == 0){
                theoryvdatarcpjetptlegend->AddEntry(FinalRCP_Pt[pT-1][cent], "STAR AuAu", "P");
                theoryvdatarcpjetptlegend->AddEntry(TheoryRCP_JetPt[pT-1][cent], "LIDO AuAu (MPI = off)", "P");
            }
            theoryvdatarcpjetptlegend->Draw("SAME");
        }
        
    }

    FinalRCP_PtCanvas_LIDO->SaveAs(Form("%s/FinalRCP_Pt_LIDO_Iter_%i.pdf", PlotDir.Data(),4));

    TCanvas *FinalRCP_ZCanvas_LIDO = new TCanvas("FinalRCP_ZCanvas_LIDO", "FinalRCP_ZCanvas_LIDO", 2100, 700);
    FinalRCP_ZCanvas_LIDO->Divide(5, 2);

    TLegend *theoryvdatarcpzlegend;
    for (int pT = 1; pT <= 5; pT++){
        if (pT != 1 && pT != 4) continue;
        for (int cent = 0; cent < 2; cent++){
            if (cent == 1) continue;
            FinalRCP_ZCanvas_LIDO->cd(cent*5+pT);
            // gPad->SetLogy();
            if (pT == 1 && cent == 0) theoryvdatarcpzlegend = new TLegend(0.1, 0.6, 0.5, 0.9);
            // FinalRCP_Z[pT-1][cent]->SetRange(5, 20);
            FinalRCP_Z_Sys[pT-1][cent]->SetLineColor(kBlue);
            FinalRCP_Z_Sys[pT-1][cent]->SetFillColor(kBlue);
            FinalRCP_Z[pT-1][cent]->SetMarkerStyle(29);
            FinalRCP_Z[pT-1][cent]->SetMarkerSize(1.2);
            FinalRCP_Z_Sys[pT-1][cent]->SetFillStyle(0);
            FinalRCP_Z_Sys[pT-1][cent]->Draw("A2");
            FinalRCP_Z[pT-1][cent]->Draw("P SAME");
            FinalRCP_Z_Sys[pT-1][cent]->GetXaxis()->SetRangeUser(0,1);
            FinalRCP_Z_Sys[pT-1][cent]->GetYaxis()->SetRangeUser(0.1, 2.5);
            TheoryRCP_JetZ[pT-1][cent]->SetLineColor(kGreen-2);
            TheoryRCP_JetZ[pT-1][cent]->SetFillColor(kGreen-2);
            TheoryRCP_JetZ[pT-1][cent]->SetMarkerColor(kGreen-2);
            TheoryRCP_JetZ[pT-1][cent]->SetMarkerStyle(20);
            TheoryRCP_JetZ[pT-1][cent]->Draw("E5 SAME");
            ZLineAtOne->Draw("SAME");
            if (pT == 1 && cent == 0){
                theoryvdatarcpzlegend->AddEntry(FinalRCP_Z[pT-1][cent], "STAR AuAu", "P");
                theoryvdatarcpzlegend->AddEntry(TheoryRCP_JetZ[pT-1][cent], "LIDO AuAu (MPI = off)", "P");
            }
            theoryvdatarcpzlegend->Draw("SAME");
        }
    }

    FinalRCP_ZCanvas_LIDO->SaveAs(Form("%s/FinalRCP_Z_LIDO_Iter_%i.pdf", PlotDir.Data(),4));

    TCanvas *FinalRCP_dRCanvas_LIDO = new TCanvas("FinalRCP_dRCanvas_LIDO", "FinalRCP_dRCanvas_LIDO", 2100, 700);
    FinalRCP_dRCanvas_LIDO->Divide(5, 2);

    TLegend *theoryvdatarcpdRlegend;
    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){
            FinalRCP_dRCanvas_LIDO->cd(cent*5+pT);
            // gPad->SetLogy();
            if (pT == 1 && cent == 0) theoryvdatarcpdRlegend = new TLegend(0.1, 0.1, 0.5, 0.4);
            // FinalRCP_dR[pT-1][cent]->SetRange(5, 20);
            FinalRCP_dR_Sys[pT-1][cent]->SetLineColor(kBlue);
            FinalRCP_dR_Sys[pT-1][cent]->SetFillColor(kBlue);
            FinalRCP_dR[pT-1][cent]->SetMarkerStyle(29);
            FinalRCP_dR[pT-1][cent]->SetMarkerSize(1.2);
            FinalRCP_dR_Sys[pT-1][cent]->SetFillStyle(0);
            FinalRCP_dR_Sys[pT-1][cent]->Draw("A2");
            FinalRCP_dR[pT-1][cent]->Draw("P SAME");
            FinalRCP_dR_Sys[pT-1][cent]->GetXaxis()->SetRangeUser(0, 0.2);
            FinalRCP_dR_Sys[pT-1][cent]->GetYaxis()->SetRangeUser(0.1, 2);
            TheoryRCP_JetdR[pT-1][cent]->SetLineColor(kGreen-2);
            TheoryRCP_JetdR[pT-1][cent]->SetFillColor(kGreen-2);
            TheoryRCP_JetdR[pT-1][cent]->SetMarkerColor(kGreen-2);
            TheoryRCP_JetdR[pT-1][cent]->SetMarkerStyle(20);
            TheoryRCP_JetdR[pT-1][cent]->Draw("E5 SAME");
            dRLineAtOne->Draw("SAME");
            if (pT == 1 && cent == 0){
                theoryvdatarcpdRlegend->AddEntry(FinalRCP_dR[pT-1][cent], "STAR AuAu", "P");
                theoryvdatarcpdRlegend->AddEntry(TheoryRCP_JetdR[pT-1][cent], "LIDO AuAu (MPI = off)", "P");
            }
            theoryvdatarcpdRlegend->Draw("SAME");
        }
    }


    TGraphAsymmErrors *gTheoryD0JetPt[5][4];
    TGraphAsymmErrors *gTheoryD0JetZ[5][4];
    TGraphAsymmErrors *gTheoryD0JetZPerJet[5][4];
    TGraphAsymmErrors *gTheoryD0JetdR[5][4];

    TGraphAsymmErrors *gTheoryRCP_JetPt[5][2];
    TGraphAsymmErrors *gTheoryRCP_JetZ[5][2];
    TGraphAsymmErrors *gTheoryRCP_JetdR[5][2];

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 4; cent++){
            // if (cent == 0 && (pT != 1 && pT != 4 && pT != 5)) {MakeTmpGraphs(gTheoryD0JetPt[pT-1][cent]); MakeTmpGraphs(gTheoryD0JetZ[pT-1][cent]); MakeTmpGraphs(gTheoryD0JetdR[pT-1][cent]); continue;}
            // if (cent == 1) {MakeTmpGraphs(gTheoryD0JetPt[pT-1][cent]); MakeTmpGraphs(gTheoryD0JetZ[pT-1][cent]); MakeTmpGraphs(gTheoryD0JetdR[pT-1][cent]); continue;}
            // if (cent == 2 && (pT != 1 && pT != 4)) {MakeTmpGraphs(gTheoryD0JetPt[pT-1][cent]); MakeTmpGraphs(gTheoryD0JetZ[pT-1][cent]); MakeTmpGraphs(gTheoryD0JetdR[pT-1][cent]); continue;}
            if (cent == 0 && (pT != 1 && pT != 4)) continue;
            if (cent == 1) continue;
            if (cent == 2 && (pT != 1 && pT != 4)) continue;
            if (cent == 3 && pT != 1) continue;
            gTheoryD0JetPt[pT-1][cent] = new TGraphAsymmErrors(TheoryD0JetPt[pT-1][cent]);
            gTheoryD0JetPt[pT-1][cent]->SetNameTitle(Form("Graph %s", TheoryD0JetPt[pT-1][cent]->GetName()), Form("Graph %s", TheoryD0JetPt[pT-1][cent]->GetName()));
            cout << TheoryD0JetPt[pT-1][cent]->GetTitle() << "\t" << gTheoryD0JetPt[pT-1][cent]->GetTitle() << endl;
            gTheoryD0JetZ[pT-1][cent] = new TGraphAsymmErrors(TheoryD0JetZ[pT-1][cent]);
            gTheoryD0JetZ[pT-1][cent]->SetNameTitle(Form("Graph %s", TheoryD0JetZ[pT-1][cent]->GetName()), Form("Graph %s", TheoryD0JetZ[pT-1][cent]->GetName()));
            cout << TheoryD0JetZ[pT-1][cent]->GetTitle() << "\t" << gTheoryD0JetZ[pT-1][cent]->GetTitle() << endl;
        }
        cout << "Here" << endl;
        for (int cent = 0; cent < 4; cent++){
            if (cent == 0 && (pT != 1 && pT != 4)) continue;
            if (cent == 1) continue;
            if (cent == 2 && (pT != 1 && pT != 4)) continue;
            if (cent == 3 && pT != 1) continue;
            gTheoryD0JetdR[pT-1][cent] = new TGraphAsymmErrors(TheoryD0JetdR[pT-1][cent]);
            gTheoryD0JetdR[pT-1][cent]->SetNameTitle(Form("Graph %s", TheoryD0JetdR[pT-1][cent]->GetName()), Form("Graph %s", TheoryD0JetdR[pT-1][cent]->GetName()));
            cout << TheoryD0JetdR[pT-1][cent]->GetTitle() << "\t" << gTheoryD0JetdR[pT-1][cent]->GetTitle() << endl;
        }
        cout << "Here" << endl;
        for (int cent = 0; cent < 2; cent++){
            // if (cent == 0 && (pT != 1 && pT != 4)) {MakeTmpGraphs(gTheoryRCP_JetPt[pT-1][cent]); continue;}
            // if (cent == 1) {MakeTmpGraphs(gTheoryRCP_JetPt[pT-1][cent]); continue;}
            if (cent == 0 && (pT != 1 && pT != 4)) continue;
            if (cent == 1) continue;
            gTheoryRCP_JetPt[pT-1][cent] = new TGraphAsymmErrors(TheoryRCP_JetPt[pT-1][cent]);
            gTheoryRCP_JetPt[pT-1][cent]->SetNameTitle(Form("Graph %s", TheoryRCP_JetPt[pT-1][cent]->GetName()), Form("Graph %s", TheoryRCP_JetPt[pT-1][cent]->GetName()));

            gTheoryRCP_JetZ[pT-1][cent] = new TGraphAsymmErrors(TheoryRCP_JetZ[pT-1][cent]);
            gTheoryRCP_JetZ[pT-1][cent]->SetNameTitle(Form("Graph %s", TheoryRCP_JetZ[pT-1][cent]->GetName()), Form("Graph %s", TheoryRCP_JetZ[pT-1][cent]->GetName()));

            gTheoryRCP_JetdR[pT-1][cent] = new TGraphAsymmErrors(TheoryRCP_JetdR[pT-1][cent]);
            gTheoryRCP_JetdR[pT-1][cent]->SetNameTitle(Form("Graph %s", TheoryRCP_JetdR[pT-1][cent]->GetName()), Form("Graph %s", TheoryRCP_JetdR[pT-1][cent]->GetName()));
        }
        cout << "Here" << endl;
    }

    TFile *FinalPlots = new TFile(Form("%s/FinalPlots.root", PlotDir.Data()), "RECREATE");
    FinalPlots->cd();
    for (int pT = 1; pT <= 5; pT++){
        FinalPlots->mkdir(Form("D0pT_%i_%i", pT, 10));
        FinalPlots->cd(Form("D0pT_%i_%i", pT, 10));
        for (int cent = 0; cent < 4; cent++){
            FinalJetPt[pT-1][cent]->Write(Form("Jet pT Spectra %s Stat", CentName[cent].Data()));
            FinalJetPt_Sys[pT-1][cent]->Write(Form("Jet pT Spectra %s Sys", CentName[cent].Data()));
            FinalJetZ[pT-1][cent]->Write(Form("Jet Z Spectra %s Stat", CentName[cent].Data()));
            FinalJetZ_Sys[pT-1][cent]->Write(Form("Jet Z Spectra %s Sys", CentName[cent].Data()));
            FinalJetZPerJet[pT-1][cent]->Write(Form("Jet Z Per Jet %s Stat", CentName[cent].Data()));
            FinalJetZPerJet_Sys[pT-1][cent]->Write(Form("Jet Z Per Jet %s Sys", CentName[cent].Data()));
            FinalJetdR[pT-1][cent]->Write(Form("Jet dR %s Stat", CentName[cent].Data()));
            FinalJetdR_Sys[pT-1][cent]->Write(Form("Jet dR %s Sys", CentName[cent].Data()));
            if (cent == 0 && (pT != 1 && pT != 4)) continue;
            if (cent == 1) continue;
            if (cent == 2 && (pT != 1 && pT != 4)) continue;
            if (cent == 3 && pT != 1) continue;
            gTheoryD0JetPt[pT-1][cent]->Write(Form("LIDO Jet pT Spectra %s Stat", CentName[cent].Data()));
            gTheoryD0JetZ[pT-1][cent]->Write(Form("LIDO Jet Z Spectra %s Stat", CentName[cent].Data()));
            gTheoryD0JetdR[pT-1][cent]->Write(Form("LIDO Jet dR %s Stat", CentName[cent].Data()));
            TheoryD0JetPt[pT-1][cent]->Write(Form("LIDO Jet pT Spectra %s Stat Hist", CentName[cent].Data()));
            TheoryD0JetZ[pT-1][cent]->Write(Form("LIDO Jet Z Spectra %s Stat Hist", CentName[cent].Data()));
            TheoryD0JetdR[pT-1][cent]->Write(Form("LIDO Jet dR %s Stat Hist", CentName[cent].Data()));
        }
        cout << "Here" << endl;
        for (int cent = 0; cent < 2; cent++){
            FinalRCP_Pt[pT-1][cent]->Write(Form("RCP Jet pT Spectra %s Stat", CentName[cent].Data()));
            FinalRCP_Pt_Sys[pT-1][cent]->Write(Form("RCP Jet pT Spectra %s Sys", CentName[cent].Data()));
            FinalRCP_Z[pT-1][cent]->Write(Form("RCP Jet Z Spectra %s Stat", CentName[cent].Data()));
            FinalRCP_Z_Sys[pT-1][cent]->Write(Form("RCP Jet Z Spectra %s Sys", CentName[cent].Data()));
            FinalRCP_ZPerJet[pT-1][cent]->Write(Form("RCP Jet Z Per Jet %s Stat", CentName[cent].Data()));
            FinalRCP_ZPerJet_Sys[pT-1][cent]->Write(Form("RCP Jet Z Per Jet %s Sys", CentName[cent].Data()));
            FinalRCP_dR[pT-1][cent]->Write(Form("RCP Jet dR %s Stat", CentName[cent].Data()));
            FinalRCP_dR_Sys[pT-1][cent]->Write(Form("RCP Jet dR %s Sys", CentName[cent].Data()));
            if (cent == 0 && (pT != 1 && pT != 4)) continue;
            if (cent == 1) continue;
            gTheoryRCP_JetPt[pT-1][cent]->Write(Form("RCP LIDO Jet pT Spectra %s Stat", CentName[cent].Data()));
            TheoryRCP_JetPt[pT-1][cent]->Write(Form("RCP LIDO Jet pT Spectra %s Stat Hist", CentName[cent].Data()));
            gTheoryRCP_JetZ[pT-1][cent]->Write(Form("RCP LIDO Jet Z Spectra %s Stat", CentName[cent].Data()));
            TheoryRCP_JetZ[pT-1][cent]->Write(Form("RCP LIDO Jet Z Spectra %s Stat Hist", CentName[cent].Data()));
            gTheoryRCP_JetdR[pT-1][cent]->Write(Form("RCP LIDO Jet dR %s Stat", CentName[cent].Data()));
            TheoryRCP_JetdR[pT-1][cent]->Write(Form("RCP LIDO Jet dR %s Stat Hist", CentName[cent].Data()));
        }
    }
    FinalPlots->Close();

    vector<TString> IterName;
    IterName.push_back("Regularisation Parameter");
    IterName.push_back("Prior Variation");
    IterName.push_back("Tracking Efficiency");
    IterName.push_back("D0 Signal Extraction");
    IterName.push_back("D0 Reconstruction No VC Efficiency");
    IterName.push_back("D0 Reconstruction VC Efficiency");

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            cout << "          ============== pT = " << pT << "-10, Cent = " << CentName[cent].Data() << " ==============" << endl;
            for (int vary = 0; vary < lowJetPt[pT-1][cent].size(); vary++){
                lowJetPt[pT-1][cent][vary]->Scale(-1.);
                cout << Form("%i %30.30s: Lower Limit = %.2f %%, Upper Limit = %.2f %%", vary+1, IterName[vary].Data(), -1.*lowJetPt[pT-1][cent][vary]->GetBinContent(lowJetPt[pT-1][cent][vary]->GetMaximumBin()), highJetPt[pT-1][cent][vary]->GetBinContent(highJetPt[pT-1][cent][vary]->GetMaximumBin())) << endl;
            }
        }
    }

    for ( int bin = 1; bin <= JetdR[0][0][0]->GetNbinsX(); bin++ ){
        cout << "Stat Errors For dR RCP" << endl;

        cout << Form("Bin # %i :::: Cent == 0, dR = %.2f - %.2f, Val = %f, Stat Err = %.2f, Stat Err Percent = %.2f", bin, JetdR[0][0][0]->GetBinLowEdge(bin), JetdR[0][0][0]->GetBinLowEdge(bin+1), JetdR[0][0][0]->GetBinContent(bin), JetdR[0][0][0]->GetBinError(bin), JetdR[0][0][0]->GetBinError(bin)/JetdR[0][0][0]->GetBinContent(bin)*100) << endl;
        // cout << Form("Cent == 1, dR = %.2f - %.2f, Stat Err = %.2f", JetdR[0][0][1]->GetBinLowEdge(bin), JetdR[0][0][1]->GetBinLowEdge(bin+1), JetdR[0][0][1]->GetBinError(bin)) << endl;
        cout << Form("Bin # %i :::: Cent == 2, dR = %.2f - %.2f, Val = %f, Stat Err = %.2f, Stat Err Percent = %.2f", bin, JetdR[0][0][2]->GetBinLowEdge(bin), JetdR[0][0][2]->GetBinLowEdge(bin+1), JetdR[0][0][2]->GetBinContent(bin), JetdR[0][0][2]->GetBinError(bin), JetdR[0][0][2]->GetBinError(bin)/JetdR[0][0][2]->GetBinContent(bin)*100) << endl;

        cout << "Stat Error For RCP From Calc" << endl;
        cout << Form("Bin # %i :::: Central = %.2f", bin, TMath::Sqrt(pow(JetdR[0][0][0]->GetBinError(bin)/JetdR[0][0][0]->GetBinContent(bin), 2) + pow(JetdR[0][0][2]->GetBinError(bin)/JetdR[0][0][2]->GetBinContent(bin), 2))) << endl;
        cout << Form("Bin # %i :::: Central From RCP Plot = %.2f", bin, RCP_dR[0][0][0]->GetBinError(bin)/RCP_dR[0][0][0]->GetBinContent(bin)) << endl;
    }
}