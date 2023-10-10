#include "Config.h"
#include <TString.h>

TFile *eff = new TFile("EssentialFiles/D0Eff.root", "READ");
TH1F *hEff0 = (TH1F *)eff->Get("hEff0");
TH1F *hEff1 = (TH1F *)eff->Get("hEff1");
TH1F *hEff2 = (TH1F *)eff->Get("hEff2");
TH1F *hEff3 = (TH1F *)eff->Get("hEff3");
TH1F *hEff4 = (TH1F *)eff->Get("hEff4");

TH1F *hEff0_Sys = (TH1F *)hEff0->Clone("hEff0_Sys");
TH1F *hEff1_Sys = (TH1F *)hEff1->Clone("hEff1_Sys");
TH1F *hEff2_Sys = (TH1F *)hEff2->Clone("hEff2_Sys");
TH1F *hEff3_Sys = (TH1F *)hEff3->Clone("hEff3_Sys");
TH1F *hEff4_Sys = (TH1F *)hEff4->Clone("hEff4_Sys");

TFile *Eff_Sys = new TFile("EssentialFiles/Sys4Matt.root", "READ");
TGraph *gSys_cent0_err0 = (TGraph *)Eff_Sys->Get("gSys_cent0_err0");
TGraph *gSys_cent0_err1 = (TGraph *)Eff_Sys->Get("gSys_cent0_err1");
TGraph *gSys_cent0_err2 = (TGraph *)Eff_Sys->Get("gSys_cent0_err2");
TGraph *gSys_cent0_err3 = (TGraph *)Eff_Sys->Get("gSys_cent0_err3");
TGraph *gSys_cent0_err4 = (TGraph *)Eff_Sys->Get("gSys_cent0_err4");
TGraph *gSys_cent0_err5 = (TGraph *)Eff_Sys->Get("gSys_cent0_err5");
TGraph *gSys_cent0_err6 = (TGraph *)Eff_Sys->Get("gSys_cent0_err6");
TGraph *gSys_cent0_err7 = (TGraph *)Eff_Sys->Get("gSys_cent0_err7");
TGraph *gSys_cent0_err8 = (TGraph *)Eff_Sys->Get("gSys_cent0_err8");
TGraph *gSys_cent0_err9 = (TGraph *)Eff_Sys->Get("gSys_cent0_err9");

TGraph *gSys_cent5_err0 = (TGraph *)Eff_Sys->Get("gSys_cent5_err0");
TGraph *gSys_cent5_err1 = (TGraph *)Eff_Sys->Get("gSys_cent5_err1");
TGraph *gSys_cent5_err2 = (TGraph *)Eff_Sys->Get("gSys_cent5_err2");
TGraph *gSys_cent5_err3 = (TGraph *)Eff_Sys->Get("gSys_cent5_err3");
TGraph *gSys_cent5_err4 = (TGraph *)Eff_Sys->Get("gSys_cent5_err4");
TGraph *gSys_cent5_err5 = (TGraph *)Eff_Sys->Get("gSys_cent5_err5");
TGraph *gSys_cent5_err6 = (TGraph *)Eff_Sys->Get("gSys_cent5_err6");
TGraph *gSys_cent5_err7 = (TGraph *)Eff_Sys->Get("gSys_cent5_err7");
TGraph *gSys_cent5_err8 = (TGraph *)Eff_Sys->Get("gSys_cent5_err8");
TGraph *gSys_cent5_err9 = (TGraph *)Eff_Sys->Get("gSys_cent5_err9");

TGraph *gSys_cent6_err0 = (TGraph *)Eff_Sys->Get("gSys_cent6_err0");
TGraph *gSys_cent6_err1 = (TGraph *)Eff_Sys->Get("gSys_cent6_err1");
TGraph *gSys_cent6_err2 = (TGraph *)Eff_Sys->Get("gSys_cent6_err2");
TGraph *gSys_cent6_err3 = (TGraph *)Eff_Sys->Get("gSys_cent6_err3");
TGraph *gSys_cent6_err4 = (TGraph *)Eff_Sys->Get("gSys_cent6_err4");
TGraph *gSys_cent6_err5 = (TGraph *)Eff_Sys->Get("gSys_cent6_err5");
TGraph *gSys_cent6_err6 = (TGraph *)Eff_Sys->Get("gSys_cent6_err6");
TGraph *gSys_cent6_err7 = (TGraph *)Eff_Sys->Get("gSys_cent6_err7");
TGraph *gSys_cent6_err8 = (TGraph *)Eff_Sys->Get("gSys_cent6_err8");
TGraph *gSys_cent6_err9 = (TGraph *)Eff_Sys->Get("gSys_cent6_err9");

TH1F *D0Sys_010 = new TH1F("D0Sys_010", "D0Sys_010", PTBINS, edges);
TH1F *D0Sys_1040 = new TH1F("D0Sys_1040", "D0Sys_1040", PTBINS, edges);
TH1F *D0Sys_4080 = new TH1F("D0Sys_480", "D0Sys_4080", PTBINS, edges);

TFile *doub = new TFile("EssentialFiles/MisPID_SB_Final.root", "READ");
TGraph *gDoub0 = (TGraph *)doub->Get("DoubleCounting_Cen_0_SB");
TGraph *gDoub1 = (TGraph *)doub->Get("DoubleCounting_Cen_1_SB");
TGraph *gDoub2 = (TGraph *)doub->Get("DoubleCounting_Cen_2_SB");
TGraph *gDoub3 = (TGraph *)doub->Get("DoubleCounting_Cen_3_SB");
bool useRefmult = true;
bool useEff = true;    // aPply eff corr if true
bool useDouble = true; // Apply double corr. if true
// bool usePt = true; // Apply Jet pT 5 3 geV if true
bool usePt = false; // Apply Jet pT > 3 geV if true
// bool useHighD0Pt = false; // Appl D0 pt > 3 GeV if true
// bool useLowD0Pt = false;// for D0 1-3 GeV //changing it to D0 1-10 GeV now

double ptlow = 1.;
double pthigh = 10.;

bool useD0Eta = false;   // D0 eta in 0.6
bool useJetArea = false; // Jet area cut
bool useBKG = false;     // FOr BKG distributions
bool useConst = false;   // (D0-Jet) Pt <0 thronw away if true
bool useZ = false;       // Cut events with D0 z <0 or > 1
bool SYS = false;        // Run systematic uncertainties
// char NAME[100] = Form("Histograms3_D0%i_%iGeV_NewBinningOldCuts.root", (int)ptlow, (int)pthigh) ;

float muCentrality;
float muSWeight;
float muRecoD0Pt;
float muRecoJetPt;
float muRecoZ;
float muRecoDeltaR;

TString NAME = Form("Histograms_D0%i_%iGeV_RecoJetPt_%i_%i.root", (int)ptlow, (int)pthigh, (int)ptlow, 1000);

void fillSys(TH1 *h, int cen, int SYSMODE)
{

    for (int i = 1; i < h->GetNbinsX() + 1; i++)
    {
        double v1, v2, v3, v4, v5, v6, v7, v8, v9, v10;
        if (cen == 1)
        {
            v1 = gSys_cent0_err0->Eval(h->GetBinCenter(i));
            v2 = gSys_cent0_err1->Eval(h->GetBinCenter(i));
            v3 = gSys_cent0_err2->Eval(h->GetBinCenter(i));
            v4 = gSys_cent0_err3->Eval(h->GetBinCenter(i));
            v5 = gSys_cent0_err4->Eval(h->GetBinCenter(i));
            v6 = gSys_cent0_err5->Eval(h->GetBinCenter(i));
            v7 = gSys_cent0_err6->Eval(h->GetBinCenter(i));
            v8 = gSys_cent0_err7->Eval(h->GetBinCenter(i));
            v9 = gSys_cent0_err8->Eval(h->GetBinCenter(i));
            v10 = gSys_cent0_err9->Eval(h->GetBinCenter(i));
        }
        else if (cen == 2)
        {
            v1 = gSys_cent5_err0->Eval(h->GetBinCenter(i));
            v2 = gSys_cent5_err1->Eval(h->GetBinCenter(i));
            v3 = gSys_cent5_err2->Eval(h->GetBinCenter(i));
            v4 = gSys_cent5_err3->Eval(h->GetBinCenter(i));
            v5 = gSys_cent5_err4->Eval(h->GetBinCenter(i));
            v6 = gSys_cent5_err5->Eval(h->GetBinCenter(i));
            v7 = gSys_cent5_err6->Eval(h->GetBinCenter(i));
            v8 = gSys_cent5_err7->Eval(h->GetBinCenter(i));
            v9 = gSys_cent5_err8->Eval(h->GetBinCenter(i));
            v10 = gSys_cent5_err9->Eval(h->GetBinCenter(i));
        }
        else if (cen == 3)
        {

            v1 = gSys_cent6_err0->Eval(h->GetBinCenter(i));
            v2 = gSys_cent6_err1->Eval(h->GetBinCenter(i));
            v3 = gSys_cent6_err2->Eval(h->GetBinCenter(i));
            v4 = gSys_cent6_err3->Eval(h->GetBinCenter(i));
            v5 = gSys_cent6_err4->Eval(h->GetBinCenter(i));
            v6 = gSys_cent6_err5->Eval(h->GetBinCenter(i));
            v7 = gSys_cent6_err6->Eval(h->GetBinCenter(i));
            v8 = gSys_cent6_err7->Eval(h->GetBinCenter(i));
            v9 = gSys_cent6_err8->Eval(h->GetBinCenter(i));
            v10 = gSys_cent6_err9->Eval(h->GetBinCenter(i));
        }
        // Ignore v3, which comes from yield extraction
        // double total = sqrt(v1*v1+v2*v2+v4*v4+v5*v5+v6*v6+v7*v7+v8*v8+v9*v9+v10*v10);
        // for ratios take out correalted ones
        // double total = sqrt(v1 * v1 + v2 * v2 + v4 * v4 + v5 * v5 + v6 * v6 + v7 * v7 + v10 * v10);
        double total = sqrt(v1 * v1 + v2 * v2 + v4 * v4 + v5 * v5 + v6 * v6 + v10 * v10); // Vertex Correction Is Looked At Separately Now. It's uncorrelated errors.
        if (abs(SYSMODE) == 1) h->SetBinContent(i, total);
        else if (abs(SYSMODE) == 2) h->SetBinContent(i,v7);
        else h->SetBinContent(i, 0);
        cout << "Cen " << cen << " bin " << i << " total " << total << endl;
    }
}

int SAVE = 1;

double getSysVarEff(int mCen, float mPt){
    double syseff = 0;
    if (mCen >= 0 && mCen < 10)
        syseff = D0Sys_010->GetBinContent(D0Sys_010->FindBin(mPt));
    else if (mCen >= 10 && mCen < 40)
        syseff = D0Sys_1040->GetBinContent(D0Sys_1040->FindBin(mPt));
    else if (mCen >= 40 && mCen < 80)
        syseff = D0Sys_4080->GetBinContent(D0Sys_4080->FindBin(mPt));
    
    return syseff;
}

double getEff(int mCen, float mPt)
{
    double eff = 0;

    // Note I am using the efficiency histograms with "_Sys" here. If running applyWeights() only, this doesn't matter
    // If running doSys(), this will shuffle the efficinecy with a Gaussian each iteration in the _Sys histograms.
    // To not have to change the code here, _Sys is used all the time but is equivelent to original histograms in nominal running
    if (mCen >= 0 && mCen < 10)
        eff = hEff0_Sys->GetBinContent(hEff0_Sys->FindBin(mPt));
    else if (mCen >= 10 && mCen < 20)
        eff = hEff1_Sys->GetBinContent(hEff1_Sys->FindBin(mPt));
    else if (mCen >= 20 && mCen < 40)
        eff = hEff2_Sys->GetBinContent(hEff2_Sys->FindBin(mPt));
    else if (mCen >= 40 && mCen < 60)
        eff = hEff3_Sys->GetBinContent(hEff3_Sys->FindBin(mPt));
    else if (mCen >= 60 && mCen < 80)
        eff = (2. / 3.) * hEff4_Sys->GetBinContent(hEff4_Sys->FindBin(mPt)); // efficiency is from paper which has scale factor
    return eff;
}
double getDoubleCount(int mCen, float mPt)
{
    double dcount = 1;
    if (mCen == 0 && mCen == 80)
        dcount = gDoub0->Eval(mPt);
    else if (mCen >= 0 && mCen < 10)
        dcount = gDoub1->Eval(mPt);
    else if (mCen >= 10 && mCen < 40)
        dcount = gDoub2->Eval(mPt);
    else if (mCen >= 40 && mCen < 80)
        dcount = gDoub3->Eval(mPt);
    return dcount;
}

void DN(TH1 *R, int col, int mar, double sum)
{
    for (int i = 1; i < R->GetNbinsX() + 1; i++)
    {
        double val = R->GetBinContent(i);
        double er = R->GetBinError(i);
        double width = R->GetBinWidth(i);
        R->SetBinContent(i, val / width / sum);
        R->SetBinError(i, er / width / sum);
    }
    //    R->Scale(1./sum);
    R->SetLineColor(col);
    R->SetMarkerColor(col);
    R->SetMarkerStyle(mar);
}


double fillHistos(char file[100], TString mode, int SYSMODE, TH1 *D0PtNoWeights, TH1 *D0Pt, TH1 *JetPt, TH1 *Z, TH1 *R, TH2 *ZPt, int cen_low, int cen_high, TH1 *hArea, TH1 *JetEta, TH1 *D0Eta, TH1 *JetPhi, TH1 *D0Phi, TH2 *D0JetEta, TH1 *Eta, TH2 *JetPt_R, TH2 *JetPt_Z, bool widebin, TH1D *JetPt_Area, TH1D *Z_Area, TH2D *ZPt_Area, TH2D *JetPt_R_Area, TH1D *JetPt_Area_Wide, TH1D *Z_Area_Wide, TH2D *ZPt_Area_Wide, TH2D *JetPt_R_Area_Wide, TH2D *JetPtAreaVsJetPtCS = NULL, TH1 *hCent = NULL, TH2 *JetAreaVsJetNConstArea = NULL, TH2 *JetAreaVsJetNConstCS = NULL, TH2 *JetAreaPtVsJetNConstArea = NULL, TH2 *JetCSPtVsJetNConstCS = NULL, TH2* JetCentVsJetNConst = NULL, TH1 *JetRefMult = NULL, TTree *outtree = NULL)
{
    if (!SYS)
        cout << " >>> On fill with: " << file << endl;
    if (outtree){
        cout << "Filling MultiFold File: " << outtree->GetName() << endl;
    }
    double sumw = 0;
    TFile *f_D = new TFile(file);
    TChain *ch = (TChain *)f_D->Get("Signal_sw");
    float mM;
    float mJA;
    float weight;
    float mPt;
    float mJPt;
    float mR;
    float mZ;
    int mCen;
    float mW;
    float mEta;
    float mJEta;
    float mPhi;
    float mJPhi;
    float mJPtArea;
    float mZArea;
    int mNConstCS;
    int mNConstArea;
    float mgRefMultCorr;
    ch->SetBranchAddress("mM", &mM);
    ch->SetBranchAddress("sWeight", &weight);
    ch->SetBranchAddress("mPt", &mPt);
    ch->SetBranchAddress("mJPt", &mJPt);
    ch->SetBranchAddress("mR", &mR);
    ch->SetBranchAddress("mZ", &mZ);
    ch->SetBranchAddress("mCen", &mCen);
    ch->SetBranchAddress("mWeight", &mW);
    ch->SetBranchAddress("mJetArea", &mJA);
    ch->SetBranchAddress("mEta", &mEta);
    ch->SetBranchAddress("mJEta", &mJEta);
    ch->SetBranchAddress("mGRefMult", &mgRefMultCorr);

    if (mode == "CS")
    {
        ch->SetBranchAddress("mJPtArea", &mJPtArea);
        ch->SetBranchAddress("mZArea", &mZArea);
        ch->SetBranchAddress("mConstCS", &mNConstCS);
        ch->SetBranchAddress("mConstArea", &mNConstArea);
    }

    double sigma = 0;
    if (cen_low == 0 && cen_high == 10)
        sigma = 20; // 0.156415;
    if (cen_low == 10 && cen_high == 20)
        sigma = 16; // 0.143552;
    if (cen_low == 20 && cen_high == 30)
        sigma = 16; // 0.135034;
    if (cen_low == 30 && cen_high == 40)
        sigma = 16; // 0.120011;
    if (cen_low == 40 && cen_high == 60)
        sigma = 10; // 0.100353;
    if (cen_low == 10 && cen_high == 40)
        sigma = 16; // 00.136703;
    if (cen_low == 40 && cen_high == 80)
        sigma = 10; // 00.0972255;
    int loop = ch->GetEntries();

    double count = 0;

    for (int i = 0; i < loop; i++)
    {
        ch->GetEntry(i);
        double eff = getEff(mCen, mPt);
        double dcount = getDoubleCount(mCen, mPt);
        double syseff = 1.;
        if (SYSMODE > 0){
            syseff = 1. + getSysVarEff(mCen, mPt);
        }
        else if (SYSMODE < 0){
            syseff = 1. - getSysVarEff(mCen, mPt);
        }

        if (mCen < cen_low || mCen >= cen_high)
            continue;
        //     if(useHighD0Pt && mPt<4)continue;
        // if(useLowD0Pt && (mPt > 5))continue;

        double mJPtAreaOld = mJPtArea;
        double mZAreaOld = mZArea;

        if (mPt < ptlow || mPt >= pthigh)
            continue;
        
        if (mode == "AreaBased"){
            if (!widebin){
                if (usePt && (mJPt < 3 || mJPt > 30))
                    continue;
            }
        }

        else{
            // if (usePt && mJPt > 30) continue;
            // mZ = mPt / mJPt;
            if (mZ >= 1.) mZ = 0.999;
            if (mJPt < mPt) mJPt = mPt;
            
            // if (mJPt > 11) mJPt = mJPtArea;
            // if (!outtree){
                
                if (mJPtArea <= nbinsjetpt_wide[0]) mJPtArea = nbinsjetpt_wide[0] + 0.001;
                if (mJPtArea >= nbinsjetpt_wide[njpt_bins_wide]) mJPtArea = nbinsjetpt_wide[njpt_bins_wide] - 0.001;
                if (mZArea <= nbinsz_wide[0]) mZArea = nbinsz_wide[0] + 0.001;
                if (mZArea >= nbinsz_wide[nz_bins_wide]) mZArea = nbinsz_wide[nz_bins_wide] - 0.001;
            // }
        }

        // if (mZ >= 1.) mZ = 0.999;
        // if (mJPt < mPt) mJPt = mPt;

        // if (mZ < 0)cout << mJPt << "\t" << mZ << "\t" << weight << endl;

        if (useD0Eta && fabs(mJEta) > 0.3)
            continue;
        if (useJetArea && mJA < 0.35)
            continue;
        if (useConst && (mPt - mJPt) > 0)
            continue;
        if (useZ && (mZ < 0. || mZ > 1.))
            continue;
        if (useEff && useDouble)
            weight *= (1. - dcount) / eff;
        else if (useDouble && !useEff)
            weight *= (1. - dcount);
        else if (!useDouble && useEff)
            weight *= 1 / eff;
        if (useRefmult)
            weight *= mW;
        if (useBKG)
            weight = 1. - weight;

        weight *= syseff;


        if (mNConstCS == 1) {
            // cout << mPt << "\t" << mJPt << "\t" << mJPtArea << "\t" << mNConstArea << endl;
            count+=weight;
        }
        // if (mCen >= 40 && mCen < 60 && mPt < 5 && mJPt >= 15 && mJPt < 38) cout << "Before = " << mCen << "\t" << mJPt << "\t" << mPt << "\t" << weight << endl;
        if (hCent) hCent->Fill(mCen, weight);
        if (D0PtNoWeights) D0PtNoWeights->Fill(mPt);
        sumw += weight;
        if (D0Pt) D0Pt->Fill(mPt, weight);
        if (JetPt) JetPt->Fill(mJPt, weight);
        if (Z) Z->Fill(mZ, weight);
        if (hArea) hArea->Fill(mJA, weight);
        if (R) R->Fill(mR, weight);
        if (ZPt) ZPt->Fill(mJPt, mZ, weight);
        if (JetEta) JetEta->Fill(mJEta, weight);
        if (D0Eta) D0Eta->Fill(mEta, weight);
        if (JetPhi) JetPhi->Fill(mJPhi, weight);
        if (D0Phi) D0Phi->Fill(mPhi, weight);
        if (D0JetEta) D0JetEta->Fill(mEta, mJEta, weight);
        if (Eta) Eta->Fill(mJEta - mEta);
        if (JetPt_R) JetPt_R->Fill(mJPt, mR, weight);
        if (JetPt_Z) JetPt_Z->Fill(mJPt, mZ, weight);

        if (mode == "CS")
        {
            if (mJPtArea >= 3 && mJPtArea <= 30){
                if (JetPt_Area) JetPt_Area->Fill(mJPtArea, weight);
                if (Z_Area) Z_Area->Fill(mZArea, weight);
                if (ZPt_Area) ZPt_Area->Fill(mJPtArea, mZArea, weight);
                if (JetPt_R_Area) JetPt_R_Area->Fill(mJPtArea, mR, weight);
            }
                
            if (JetPt_Area_Wide) JetPt_Area_Wide->Fill(mJPtArea, weight);
            if (Z_Area_Wide) Z_Area_Wide->Fill(mZArea, weight);
            if (ZPt_Area_Wide) ZPt_Area_Wide->Fill(mJPtArea, mZArea, weight);
            // cout << JetPt_R_Area_Wide->GetName() << "\t" << JetPt_R_Area_Wide->Fill(mJPtArea, mR, weight) << endl;
            if (JetPt_R_Area_Wide) JetPt_R_Area_Wide->Fill(mJPtArea, mR, weight);

            if (JetPtAreaVsJetPtCS) JetPtAreaVsJetPtCS->Fill(mJPtArea, mJPt, weight);
            if (JetAreaVsJetNConstArea) JetAreaVsJetNConstArea->Fill(mJA, mNConstArea, weight);
            if (JetAreaVsJetNConstCS) JetAreaVsJetNConstCS->Fill(mJA, mNConstCS, weight);
            if (JetAreaPtVsJetNConstArea) JetAreaPtVsJetNConstArea->Fill(mJPtArea, mNConstArea, weight);
            if (JetCSPtVsJetNConstCS) JetCSPtVsJetNConstCS->Fill(mJPt, mNConstCS, weight);
            if (JetCentVsJetNConst) JetCentVsJetNConst->Fill(mCen, mNConstCS, weight);
            if (JetRefMult) JetRefMult->Fill(mgRefMultCorr, weight);

            if (outtree){
                // cout << "Gets here" << endl;
                muCentrality = mCen;
                muSWeight = weight;
                muRecoD0Pt = mPt;
                muRecoJetPt = mJPtAreaOld;
                muRecoZ = mZAreaOld;
                muRecoDeltaR = mR;
                outtree->Fill();
            }

            // cout << mJPtArea << "\t" << mJPt << "\t" << weight << endl;
        }

        // if (mZ < 0)cout << mJPt << "\t" << mZ << "\t" << weight << endl;
    }
    f_D->Close();
    // cout << "Integral From Step = " << sumw << "\t" << ZPt->Integral() << "\t" << JetPt_Area_Wide->Integral() << "\t" << JetPt_R_Area_Wide->Integral() << endl;
    
    // cout << "Single Const. Jet Ratio = " << count << "\t" << sumw << "\t" << count/sumw << endl;
    return sumw;
}

void shuffleEff()
{
    for (int i = 1; i < hEff0->GetNbinsX() + 1; i++)
    {
        double val0 = hEff0->GetBinContent(i);
        double val1 = hEff1->GetBinContent(i);
        double val2 = hEff2->GetBinContent(i);
        double val3 = hEff3->GetBinContent(i);
        double val4 = hEff4->GetBinContent(i);
        double err0 = D0Sys_010->GetBinContent(i);  // hEff0->GetBinError(i);
        double err1 = D0Sys_1040->GetBinContent(i); // hEff1->GetBinError(i);
        double err2 = D0Sys_1040->GetBinContent(i); // hEff2->GetBinError(i);
        double err3 = D0Sys_1040->GetBinContent(i); // hEff3->GetBinError(i);
        double err4 = D0Sys_4080->GetBinContent(i); // hEff4->GetBinError(i);
        float shift0 = gRandom->Gaus(0, err0);
        float shift1 = gRandom->Gaus(0, err1);
        float shift2 = gRandom->Gaus(0, err2);
        float shift3 = gRandom->Gaus(0, err3);
        float shift4 = gRandom->Gaus(0, err4);
        hEff0_Sys->SetBinContent(i, val0 / (1 + shift0));
        hEff1_Sys->SetBinContent(i, val1 / (1 + shift1));
        hEff2_Sys->SetBinContent(i, val2 / (1 + shift2));
        hEff3_Sys->SetBinContent(i, val3 / (1 + shift3));
        hEff4_Sys->SetBinContent(i, val4 / (1 + shift4));
    }
}

void RemoveNegativeYieldBins(TH1 *h)
{
    for (int i = 0; i <= h->GetNbinsX() + 1; i++)
    {
        if (h->GetBinContent(i) < 0)
        {
            h->SetBinContent(i, 0);
            h->SetBinError(i, h->GetBinError(i));
        }
    }
}

void RemoveNegativeYieldBins(TH2 *h)
{
    for (int i = 0; i < h->GetNbinsX() + 1; i++)
    {
        for (int j = 0; j < h->GetNbinsY() + 1; j++)
        {
            if (h->GetBinContent(i, j) < 0)
            {
                h->SetBinContent(i, j, 0);
                h->SetBinError(i, j, h->GetBinError(i, j));
            }
        }
    }
}

void applyweights(int SYSMODE = 0, int save = SAVE, TString mode = "CS", bool widebin = true)
{
    if (mode == "CS") widebin = false;
    // if (!SYS)
    //     gROOT->ProcessLine(".x ../myStyle.C");

    fillSys(D0Sys_010, 1, SYSMODE);
    fillSys(D0Sys_1040, 2, SYSMODE);
    fillSys(D0Sys_4080, 3, SYSMODE);
    
    if (SYSMODE == 1){
        NAME.ReplaceAll(".root", "_NoVCSysUp.root");
    }
    else if (SYSMODE == -1){
        NAME.ReplaceAll(".root", "_NoVCSysDown.root");
    }
    else if (SYSMODE == 2){
        NAME.ReplaceAll(".root", "_VCSysUp.root");
    }
    else if (SYSMODE == -2){
        NAME.ReplaceAll(".root", "_VCSysDown.root");
    }

    cout << "===================="
         << "D0 pT Bin = " << ptlow << "\t" << pthigh << "====================" << endl;
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    double sumw = 0;
    double sumw1 = 0;
    double sumw2 = 0;
    double sumw3 = 0;
    double sumw4 = 0;
    double sumw5 = 0;
    double sumw6 = 0;
    double sumw7 = 0;
    double binning[6] = {0, 0.05, 0.1, 0.2, 0.3, 0.4};
    double binning1[9] = {5, 7, 9, 11, 13, 15, 20, 30, 40};
    double binning2[11] = {0, 1, 2, 3, 4, 5, 6, 8, 10};
    TH2F *PT2D = new TH2F("PT2D", "PT2D", 35, 5, 40, 50, 0, 10); // 8,binning1,10,binning2);
    TH1F *hArea = new TH1F("hArea", "hArea", 100, 0, 1);
    TH2F *D0JPt_5 = new TH2F("D0JPt_5", "D0JPt_5", nbins_jpt, binning_jpt, PTBINS, edges); // 50,0,50);
    TH2F *Hard2D2_5 = new TH2F("r_Hard2D2_5", "r_Hard2D2_5", 100, 0, 50, 100, 0, 5);
    
    TH1F *Eta = new TH1F("Eta", "D0Eta", 200, -1, 1);
    TH1F *Eta_1 = new TH1F("Eta_1", "Eta_1", 200, -1, 1);
    TH1F *Eta_2 = new TH1F("Eta_2", "Eta_2", 200, -1, 1);
    TH1F *Eta_3 = new TH1F("Eta_3", "Eta_3", 200, -1, 1);
    TH1F *Eta_4 = new TH1F("Eta_4", "Eta_4", 200, -1, 1);
    TH1F *Eta_5 = new TH1F("Eta_5", "Eta_5", 200, -1, 1);
    TH1F *Eta_6 = new TH1F("Eta_6", "Eta_6", 200, -1, 1);
    TH1F *Eta_7 = new TH1F("Eta_7", "Eta_7", 200, -1, 1);

    TH1F *D0Eta = new TH1F("D0Eta", "D0Eta", 200, -1, 1);
    TH1F *D0Eta_1 = new TH1F("D0Eta_1", "D0Eta_1", 200, -1, 1);
    TH1F *D0Eta_2 = new TH1F("D0Eta_2", "D0Eta_2", 200, -1, 1);
    TH1F *D0Eta_3 = new TH1F("D0Eta_3", "D0Eta_3", 200, -1, 1);
    TH1F *D0Eta_4 = new TH1F("D0Eta_4", "D0Eta_4", 200, -1, 1);
    TH1F *D0Eta_5 = new TH1F("D0Eta_5", "D0Eta_5", 200, -1, 1);
    TH1F *D0Eta_6 = new TH1F("D0Eta_6", "D0Eta_6", 200, -1, 1);
    TH1F *D0Eta_7 = new TH1F("D0Eta_7", "D0Eta_7", 200, -1, 1);

    TH1F *JetEta = new TH1F("JetEta", "JetEta", 200, -1, 1);
    TH1F *JetEta_1 = new TH1F("JetEta_1", "JetEta_1", 200, -1, 1);
    TH1F *JetEta_2 = new TH1F("JetEta_2", "JetEta_2", 200, -1, 1);
    TH1F *JetEta_3 = new TH1F("JetEta_3", "JetEta_3", 200, -1, 1);
    TH1F *JetEta_4 = new TH1F("JetEta_4", "JetEta_4", 200, -1, 1);
    TH1F *JetEta_5 = new TH1F("JetEta_5", "JetEta_5", 200, -1, 1);
    TH1F *JetEta_6 = new TH1F("JetEta_6", "JetEta_6", 200, -1, 1);
    TH1F *JetEta_7 = new TH1F("JetEta_7", "JetEta_7", 200, -1, 1);

    TH2F *D0JetEta = new TH2F("D0JetEta", "D0JetEta", 200, -1, 1, 200, -1, 1);
    TH2F *D0JetEta_1 = new TH2F("D0JetEta_1", "D0JetEta_1", 200, -1, 1, 200, -1, 1);
    TH2F *D0JetEta_2 = new TH2F("D0JetEta_2", "D0JetEta_2", 200, -1, 1, 200, -1, 1);
    TH2F *D0JetEta_3 = new TH2F("D0JetEta_3", "D0JetEta_3", 200, -1, 1, 200, -1, 1);
    TH2F *D0JetEta_4 = new TH2F("D0JetEta_4", "D0JetEta_4", 200, -1, 1, 200, -1, 1);
    TH2F *D0JetEta_5 = new TH2F("D0JetEta_5", "D0JetEta_5", 200, -1, 1, 200, -1, 1);
    TH2F *D0JetEta_6 = new TH2F("D0JetEta_6", "D0JetEta_6", 200, -1, 1, 200, -1, 1);
    TH2F *D0JetEta_7 = new TH2F("D0JetEta_7", "D0JetEta_7", 200, -1, 1, 200, -1, 1);

    TH1F *D0Phi = new TH1F("D0Phi", "D0Phi", 200, -1, 1);
    TH1F *D0Phi_1 = new TH1F("D0Phi_1", "D0Phi_1", 200, -1, 1);
    TH1F *D0Phi_2 = new TH1F("D0Phi_2", "D0Phi_2", 200, -1, 1);
    TH1F *D0Phi_3 = new TH1F("D0Phi_3", "D0Phi_3", 200, -1, 1);
    TH1F *D0Phi_4 = new TH1F("D0Phi_4", "D0Phi_4", 200, -1, 1);
    TH1F *D0Phi_5 = new TH1F("D0Phi_5", "D0Phi_5", 200, -1, 1);
    TH1F *D0Phi_6 = new TH1F("D0Phi_6", "D0Phi_6", 200, -1, 1);
    TH1F *D0Phi_7 = new TH1F("D0Phi_7", "D0Phi_7", 200, -1, 1);

    TH1F *JetPhi = new TH1F("JetPhi", "JetPhi", 200, -1, 1);
    TH1F *JetPhi_1 = new TH1F("JetPhi_1", "JetPhi_1", 200, -1, 1);
    TH1F *JetPhi_2 = new TH1F("JetPhi_2", "JetPhi_2", 200, -1, 1);
    TH1F *JetPhi_3 = new TH1F("JetPhi_3", "JetPhi_3", 200, -1, 1);
    TH1F *JetPhi_4 = new TH1F("JetPhi_4", "JetPhi_4", 200, -1, 1);
    TH1F *JetPhi_5 = new TH1F("JetPhi_5", "JetPhi_5", 200, -1, 1);
    TH1F *JetPhi_6 = new TH1F("JetPhi_6", "JetPhi_6", 200, -1, 1);
    TH1F *JetPhi_7 = new TH1F("JetPhi_7", "JetPhi_7", 200, -1, 1);

    // TH1F *hArea = new TH1F("hArea", "hArea", 100, 0, 1);
    TH1F *hArea_1 = new TH1F("hArea_1", "hArea_1", 100, 0, 1);
    TH1F *hArea_2 = new TH1F("hArea_2", "hArea_2", 100, 0, 1);
    TH1F *hArea_3 = new TH1F("hArea_3", "hArea_3", 100, 0, 1);
    TH1F *hArea_4 = new TH1F("hArea_4", "hArea_4", 100, 0, 1);
    TH1F *hArea_5 = new TH1F("hArea_5", "hArea_5", 100, 0, 1);
    TH1F *hArea_6 = new TH1F("hArea_6", "hArea_6", 100, 0, 1);
    TH1F *hArea_7 = new TH1F("hArea_7", "hArea_7", 100, 0, 1);

    TH1F *centrality = new TH1F("centrality", "centrality", 20, 0, 100);

    const int nd0ptbins = 5;
    int D0PtBin[nd0ptbins + 1] = {1, 2, 3, 4, 5, 10};

    // 2D Array of 8 x ND0PtBins histograms for D0PtUnweighted, D0Pt, JetPt, Z, R, ZPt, JetPt_R
    TH1D *D0PtUnweighted[8][nd0ptbins];
    TH1D *D0Pt[8][nd0ptbins];
    TH1D *R[8][nd0ptbins];
    TH1D *JetPt[8][nd0ptbins];
    TH1D *Z[8][nd0ptbins];
    TH2D *ZPt[8][nd0ptbins];
    TH2D *JetPt_R[8][nd0ptbins];
    TH2D *JetPt_Z[8][nd0ptbins];

    TH1D *JetPt_Area[8][nd0ptbins];
    TH1D *Z_Area[8][nd0ptbins];
    TH2D *ZPt_Area[8][nd0ptbins];
    TH2D *JetPt_R_Area[8][nd0ptbins];

    TH1D *JetPt_Area_Wide[8][nd0ptbins];
    TH1D *Z_Area_Wide[8][nd0ptbins];
    TH2D *ZPt_Area_Wide[8][nd0ptbins];
    TH2D *JetPt_R_Area_Wide[8][nd0ptbins];

    TH2D *JetPtAreaVsJetPtCS[8][nd0ptbins];

    TH2D *JetAreaVsJetNConstArea[8][nd0ptbins];
    TH2D *JetAreaVsJetNConstCS[8][nd0ptbins];

    TH2D *JetPtAreaVsJetNConstArea[8][nd0ptbins];
    TH2D *JetPtCSVsJetNConstCS[8][nd0ptbins];

    TH2D *JetCentVsJetNConst[8][nd0ptbins];

    TH1D *JetRefMult[8][nd0ptbins];

    // Define the histograms declared above based on the old definitions
    for (int i = 0; i < 8; i++)
    { // Centrality Bins
        for (int bin = 0; bin < nd0ptbins; bin++)
        { // D0 pT Bins
            D0PtUnweighted[i][bin] = new TH1D(Form("D0PtUnweighted_%i_%i", i, bin), Form("D0PtUnweighted_%i_%i", i, bin), PTBINS, edges);
            D0Pt[i][bin] = new TH1D(Form("D0Pt_%i_%i", i, bin), Form("D0Pt_%i_%i", i, bin), PTBINS, edges);
            R[i][bin] = new TH1D(Form("R_%i_%i", i, bin), Form("R_%i_%i", i, bin), nBinsdR, dRBins);
            JetPt_Z[i][bin] = new TH2D(Form("JetPt_Z_%i_%i", i, bin), Form("JetPt_Z_%i_%i", i, bin), 71, -30.5, 40.5, 2000, -1000, 1000);

            if (mode == "AreaBased"){
                if (!widebin){
                    JetPt[i][bin] = new TH1D(Form("JetPt_%i_%i", i, bin), Form("JetPt_%i_%i", i, bin), njpt_bins, nbinsjetpt);
                    Z[i][bin] = new TH1D(Form("Z_%i_%i", i, bin), Form("Z_%i_%i", i, bin), nz_bins, nbinsz);
                    ZPt[i][bin] = new TH2D(Form("ZPt_%i_%i", i, bin), Form("ZPt_%i_%i", i, bin), njpt_bins, nbinsjetpt, nz_bins, nbinsz);
                    JetPt_R[i][bin] = new TH2D(Form("JetPt_R_%i_%i", i, bin), Form("JetPt_R_%i_%i", i, bin), njpt_bins, nbinsjetpt, nBinsdR, dRBins);          
                }
                else{
                    JetPt[i][bin] = new TH1D(Form("JetPt_%i_%i", i, bin), Form("JetPt_%i_%i", i, bin), njpt_bins_wide, nbinsjetpt_wide);
                    Z[i][bin] = new TH1D(Form("Z_%i_%i", i, bin), Form("Z_%i_%i", i, bin), nz_bins_wide, nbinsz_wide);
                    ZPt[i][bin] = new TH2D(Form("ZPt_%i_%i", i, bin), Form("ZPt_%i_%i", i, bin), njpt_bins_wide, nbinsjetpt_wide, nz_bins_wide, nbinsz_wide);
                    JetPt_R[i][bin] = new TH2D(Form("JetPt_R_%i_%i", i, bin), Form("JetPt_R_%i_%i", i, bin), njpt_bins_wide, nbinsjetpt_wide, nBinsdR, dRBins);          
                }

            }
            else if (mode == "CS"){
                JetPt[i][bin] = new TH1D(Form("JetPt_%i_%i", i, bin), Form("JetPt_%i_%i", i, bin), njpt_bins_cs, nbinsjetptcs);
                Z[i][bin] = new TH1D(Form("Z_%i_%i", i, bin), Form("Z_%i_%i", i, bin), nz_bins_cs, nbinszcs);
                ZPt[i][bin] = new TH2D(Form("ZPt_%i_%i", i, bin), Form("ZPt_%i_%i", i, bin), njpt_bins_cs, nbinsjetptcs, nz_bins_cs, nbinszcs);
                JetPt_R[i][bin] = new TH2D(Form("JetPt_R_%i_%i", i, bin), Form("JetPt_R_%i_%i", i, bin), njpt_bins_cs, nbinsjetptcs, nBinsdR, dRBins);

                JetPt_Area[i][bin] = new TH1D(Form("JetPt_Area_%i_%i", i, bin), Form("JetPt_Area_%i_%i", i, bin), njpt_bins, nbinsjetpt);
                Z_Area[i][bin] = new TH1D(Form("Z_Area_%i_%i", i, bin), Form("Z_Area_%i_%i", i, bin), nz_bins, nbinsz);
                ZPt_Area[i][bin] = new TH2D(Form("ZPt_Area_%i_%i", i, bin), Form("ZPt_Area_%i_%i", i, bin), njpt_bins, nbinsjetpt, nz_bins, nbinsz);
                JetPt_R_Area[i][bin] = new TH2D(Form("JetPt_R_Area_%i_%i", i, bin), Form("JetPt_R_Area_%i_%i", i, bin), njpt_bins, nbinsjetpt, nBinsdR, dRBins);

                JetPt_Area_Wide[i][bin] = new TH1D(Form("JetPt_Area_Wide_%i_%i", i, bin), Form("JetPt_Area_Wide_%i_%i", i, bin), njpt_bins_wide, nbinsjetpt_wide);
                Z_Area_Wide[i][bin] = new TH1D(Form("Z_Area_Wide_%i_%i", i, bin), Form("Z_Area_Wide_%i_%i", i, bin), nz_bins_wide, nbinsz_wide);
                ZPt_Area_Wide[i][bin] = new TH2D(Form("ZPt_Area_Wide_%i_%i", i, bin), Form("ZPt_Area_Wide_%i_%i", i, bin), njpt_bins_wide, nbinsjetpt_wide, nz_bins_wide, nbinsz_wide);
                JetPt_R_Area_Wide[i][bin] = new TH2D(Form("JetPt_R_Area_Wide_%i_%i", i, bin), Form("JetPt_R_Area_Wide_%i_%i", i, bin), njpt_bins_wide, nbinsjetpt_wide, nBinsdR, dRBins);

                JetPtAreaVsJetPtCS[i][bin] = new TH2D(Form("JetPtAreaVsJetPtCS_%i_%i", i, bin), Form("JetPtAreaVsJetPtCS_%i_%i", i, bin), 38, -25.5, 50.5, 38, -25.5, 50.5); 
                JetAreaVsJetNConstArea[i][bin] = new TH2D(Form("JetAreaVsJetNConstArea_%i_%i", i, bin), Form("JetAreaVsJetNConstArea_%i_%i", i, bin), 100, 0., 1., 151, -0.5, 150.5);
                JetAreaVsJetNConstCS[i][bin] = new TH2D(Form("JetAreaVsJetNConstCS_%i_%i", i, bin), Form("JetAreaVsJetNConstCS_%i_%i", i, bin), 100, 0., 1., 151, -0.5, 150.5);
                
                JetPtAreaVsJetNConstArea[i][bin] = new TH2D(Form("JetPtAreaVsJetNConstArea_%i_%i", i, bin), Form("JetPtAreaVsJetNConstArea_%i_%i", i, bin), njpt_bins_wide, nbinsjetpt_wide, 151, -0.5, 150.5);
                JetPtCSVsJetNConstCS[i][bin] = new TH2D(Form("JetPtCSVsJetNConstCS_%i_%i", i, bin), Form("JetPtCSVsJetNConstCS_%i_%i", i, bin), njpt_bins_cs, nbinsjetptcs, 151, -0.5, 150.5);

                JetCentVsJetNConst[i][bin] = new TH2D(Form("JetCentVsJetNConst_%i_%i", i, bin), Form("JetCentVsJetNConst_%i_%i", i, bin), nCentBins, CentBins, njetconstbin, jetconstbin);

                JetRefMult[i][bin] = new TH1D(Form("JetRefMult_%i_%i", i, bin), Form("JetRefMult_%i_%i", i, bin), 100, 0, 800);
            }
        }
    } 

    TFile *outformultifold = new TFile(Form("MultifoldInput_D0Pt_%i_%i.root", (int)ptlow, (int)pthigh), "RECREATE");
    TTree *MuTreeForMultiFold = new TTree("MuTree", "MuTree");
    MuTreeForMultiFold->SetDirectory(outformultifold);

    MuTreeForMultiFold->Branch("Centrality", &muCentrality, "muCentrality/F");
    MuTreeForMultiFold->Branch("SWeight", &muSWeight, "muSWeight/F");
    MuTreeForMultiFold->Branch("RecoD0Pt", &muRecoD0Pt, "muRecoD0Pt/F");
    MuTreeForMultiFold->Branch("RecoJetPt", &muRecoJetPt, "muRecoJetPt/F");
    MuTreeForMultiFold->Branch("RecoZ", &muRecoZ, "muRecoZ/F");
    MuTreeForMultiFold->Branch("RecoDeltaR", &muRecoDeltaR, "muRecoDeltaR/F"); 

    // for (int bin = 0; bin < 1; bin++)
    // {
    //     // // The ones above this line are for QA only.
    //     // sumw += fillHistos(Form("./sWeights/NewsWeight_PT_%i_%i_Cen_0_80.root", 1, 5), mode, D0PtUnweighted[0][bin], D0Pt[0][bin], JetPt[0][bin], Z[0][bin], R[0][bin], ZPt[0][bin], 0, 80, hArea, JetEta, D0Eta, JetPhi, D0Phi, D0JetEta, Eta, JetPt_R[0][bin], JetPt_Z[0][bin], widebin, JetPt_Area[0][bin], Z_Area[0][bin], ZPt_Area[0][bin], JetPt_R_Area[0][bin], JetPt_Area_Wide[0][bin], Z_Area_Wide[0][bin], ZPt_Area_Wide[0][bin], JetPt_R_Area_Wide[0][bin], JetPtAreaVsJetPtCS[0][bin], centrality, JetAreaVsJetNConstArea[0][bin], JetAreaVsJetNConstCS[0][bin], JetPtAreaVsJetNConstArea[0][bin], JetPtCSVsJetNConstCS[0][bin], JetCentVsJetNConst[0][bin], JetRefMult[0][bin]);
    //     // sumw1 += fillHistos(Form("./sWeights/NewsWeight_PT_%i_%i_Cen_0_10.root", 1, 5), mode, D0PtUnweighted[1][bin], D0Pt[1][bin], JetPt[1][bin], Z[1][bin], R[1][bin], ZPt[1][bin], 0, 10, hArea_1, JetEta_1, D0Eta_1, JetPhi_1, D0Phi_1, D0JetEta_1, Eta_1, JetPt_R[1][bin], JetPt_Z[1][bin], widebin, JetPt_Area[1][bin], Z_Area[1][bin], ZPt_Area[1][bin], JetPt_R_Area[1][bin], JetPt_Area_Wide[1][bin], Z_Area_Wide[1][bin], ZPt_Area_Wide[1][bin], JetPt_R_Area_Wide[1][bin], JetPtAreaVsJetPtCS[1][bin], NULL, JetAreaVsJetNConstArea[1][bin], JetAreaVsJetNConstCS[1][bin], JetPtAreaVsJetNConstArea[1][bin], JetPtCSVsJetNConstCS[1][bin]);
    //     // sumw2 += fillHistos(Form("./sWeights/NewsWeight_PT_%i_%i_Cen_10_20.root", 1, 5), mode, D0PtUnweighted[2][bin], D0Pt[2][bin], JetPt[2][bin], Z[2][bin], R[2][bin], ZPt[2][bin], 10, 20, hArea_2, JetEta_2, D0Eta_2, JetPhi_2, D0Phi_2, D0JetEta_2, Eta_2, JetPt_R[2][bin], JetPt_Z[2][bin], widebin, JetPt_Area[2][bin], Z_Area[2][bin], ZPt_Area[2][bin], JetPt_R_Area[2][bin], JetPt_Area_Wide[2][bin], Z_Area_Wide[2][bin], ZPt_Area_Wide[2][bin], JetPt_R_Area_Wide[2][bin], JetPtAreaVsJetPtCS[2][bin], NULL, JetAreaVsJetNConstArea[2][bin], JetAreaVsJetNConstCS[2][bin], JetPtAreaVsJetNConstArea[2][bin], JetPtCSVsJetNConstCS[2][bin]);
    //     // sumw3 += fillHistos(Form("./sWeights/NewsWeight_PT_%i_%i_Cen_20_30.root", 1, 5), mode, D0PtUnweighted[3][bin], D0Pt[3][bin], JetPt[3][bin], Z[3][bin], R[3][bin], ZPt[3][bin], 20, 30, hArea_3, JetEta_3, D0Eta_3, JetPhi_3, D0Phi_3, D0JetEta_3, Eta_3, JetPt_R[3][bin], JetPt_Z[3][bin], widebin, JetPt_Area[3][bin], Z_Area[3][bin], ZPt_Area[3][bin], JetPt_R_Area[3][bin], JetPt_Area_Wide[3][bin], Z_Area_Wide[3][bin], ZPt_Area_Wide[3][bin], JetPt_R_Area_Wide[3][bin], JetPtAreaVsJetPtCS[3][bin], NULL, JetAreaVsJetNConstArea[3][bin], JetAreaVsJetNConstCS[3][bin], JetPtAreaVsJetNConstArea[3][bin], JetPtCSVsJetNConstCS[3][bin]);
    //     // sumw4 += fillHistos(Form("./sWeights/NewsWeight_PT_%i_%i_Cen_30_40.root", 1, 5), mode, D0PtUnweighted[4][bin], D0Pt[4][bin], JetPt[4][bin], Z[4][bin], R[4][bin], ZPt[4][bin], 30, 40, hArea_4, JetEta_4, D0Eta_4, JetPhi_4, D0Phi_4, D0JetEta_4, Eta_4, JetPt_R[4][bin], JetPt_Z[4][bin], widebin, JetPt_Area[4][bin], Z_Area[4][bin], ZPt_Area[4][bin], JetPt_R_Area[4][bin], JetPt_Area_Wide[4][bin], Z_Area_Wide[4][bin], ZPt_Area_Wide[4][bin], JetPt_R_Area_Wide[4][bin], JetPtAreaVsJetPtCS[4][bin], NULL, JetAreaVsJetNConstArea[4][bin], JetAreaVsJetNConstCS[4][bin], JetPtAreaVsJetNConstArea[4][bin], JetPtCSVsJetNConstCS[4][bin]);
    //     // sumw6 += fillHistos(Form("./sWeights/NewsWeight_PT_%i_%i_Cen_10_20.root", 1, 5), mode, D0PtUnweighted[6][bin], D0Pt[6][bin], JetPt[6][bin], Z[6][bin], R[6][bin], ZPt[6][bin], 10, 20, hArea_6, JetEta_6, D0Eta_6, JetPhi_6, D0Phi_6, D0JetEta_6, Eta_6, JetPt_R[6][bin], JetPt_Z[6][bin], widebin, JetPt_Area[6][bin], Z_Area[6][bin], ZPt_Area[6][bin], JetPt_R_Area[6][bin], JetPt_Area_Wide[6][bin], Z_Area_Wide[6][bin], ZPt_Area_Wide[6][bin], JetPt_R_Area_Wide[6][bin], JetPtAreaVsJetPtCS[6][bin], NULL, JetAreaVsJetNConstArea[6][bin], JetAreaVsJetNConstCS[6][bin], JetPtAreaVsJetNConstArea[6][bin], JetPtCSVsJetNConstCS[6][bin]);
    //     // sumw6 += fillHistos(Form("./sWeights/NewsWeight_PT_%i_%i_Cen_20_30.root", 1, 5), mode, D0PtUnweighted[6][bin], D0Pt[6][bin], JetPt[6][bin], Z[6][bin], R[6][bin], ZPt[6][bin], 20, 30, hArea_6, JetEta_6, D0Eta_6, JetPhi_6, D0Phi_6, D0JetEta_6, Eta_6, JetPt_R[6][bin], JetPt_Z[6][bin], widebin, JetPt_Area[6][bin], Z_Area[6][bin], ZPt_Area[6][bin], JetPt_R_Area[6][bin], JetPt_Area_Wide[6][bin], Z_Area_Wide[6][bin], ZPt_Area_Wide[6][bin], JetPt_R_Area_Wide[6][bin], JetPtAreaVsJetPtCS[6][bin], NULL, JetAreaVsJetNConstArea[6][bin], JetAreaVsJetNConstCS[6][bin], JetPtAreaVsJetNConstArea[6][bin], JetPtCSVsJetNConstCS[6][bin]);
    //     // sumw6 += fillHistos(Form("./sWeights/NewsWeight_PT_%i_%i_Cen_30_40.root", 1, 5), mode, D0PtUnweighted[6][bin], D0Pt[6][bin], JetPt[6][bin], Z[6][bin], R[6][bin], ZPt[6][bin], 30, 40, hArea_6, JetEta_6, D0Eta_6, JetPhi_6, D0Phi_6, D0JetEta_6, Eta_6, JetPt_R[6][bin], JetPt_Z[6][bin], widebin, JetPt_Area[6][bin], Z_Area[6][bin], ZPt_Area[6][bin], JetPt_R_Area[6][bin], JetPt_Area_Wide[6][bin], Z_Area_Wide[6][bin], ZPt_Area_Wide[6][bin], JetPt_R_Area_Wide[6][bin], JetPtAreaVsJetPtCS[6][bin], NULL, JetAreaVsJetNConstArea[6][bin], JetAreaVsJetNConstCS[6][bin], JetPtAreaVsJetNConstArea[6][bin], JetPtCSVsJetNConstCS[6][bin]);
        
    //     // sumw7 += fillHistos(Form("./sWeights/NewsWeight_PT_%i_%i_Cen_40_80.root", 1, 5), mode, D0PtUnweighted[7][bin], D0Pt[7][bin], JetPt[7][bin], Z[7][bin], R[7][bin], ZPt[7][bin], 40, 80, hArea_7, JetEta_7, D0Eta_7, JetPhi_7, D0Phi_7, D0JetEta_7, Eta_7, JetPt_R[7][bin], JetPt_Z[7][bin], widebin, JetPt_Area[7][bin], Z_Area[7][bin], ZPt_Area[7][bin], JetPt_R_Area[7][bin], JetPt_Area_Wide[7][bin], Z_Area_Wide[7][bin], ZPt_Area_Wide[7][bin], JetPt_R_Area_Wide[7][bin], JetPtAreaVsJetPtCS[7][bin], NULL, JetAreaVsJetNConstArea[7][bin], JetAreaVsJetNConstCS[7][bin], JetPtAreaVsJetNConstArea[7][bin], JetPtCSVsJetNConstCS[7][bin]);

    //     // cout << "Sum = " << sumw1 << "\t" << sumw6 << "\t" << sumw7 << endl;
    //     // cout << "Sum From hist = " << D0Pt[1][0]->Integral() << "\t" << D0Pt[6][0]->Integral() << "\t" << D0Pt[7][0]->Integral() << endl;

    //     sumw += fillHistos(Form("./sWeights/NewsWeight_PT_%i_%i_Cen_0_80.root", (int)ptlow, (int)pthigh), mode, D0PtUnweighted[0][bin], D0Pt[0][bin], JetPt[0][bin], Z[0][bin], R[0][bin], ZPt[0][bin], 0, 80, hArea, JetEta, D0Eta, JetPhi, D0Phi, D0JetEta, Eta, JetPt_R[0][bin], JetPt_Z[0][bin], widebin, JetPt_Area[0][bin], Z_Area[0][bin], ZPt_Area[0][bin], JetPt_R_Area[0][bin], JetPt_Area_Wide[0][bin], Z_Area_Wide[0][bin], ZPt_Area_Wide[0][bin], JetPt_R_Area_Wide[0][bin], JetPtAreaVsJetPtCS[0][bin], centrality, JetAreaVsJetNConstArea[0][bin], JetAreaVsJetNConstCS[0][bin], JetPtAreaVsJetNConstArea[0][bin], JetPtCSVsJetNConstCS[0][bin], JetCentVsJetNConst[0][bin], JetRefMult[0][bin]);
    //     sumw1 += fillHistos(Form("./sWeights/NewsWeight_PT_%i_%i_Cen_0_10.root", (int)ptlow, (int)pthigh), mode, D0PtUnweighted[1][bin], D0Pt[1][bin], JetPt[1][bin], Z[1][bin], R[1][bin], ZPt[1][bin], 0, 10, hArea_1, JetEta_1, D0Eta_1, JetPhi_1, D0Phi_1, D0JetEta_1, Eta_1, JetPt_R[1][bin], JetPt_Z[1][bin], widebin, JetPt_Area[1][bin], Z_Area[1][bin], ZPt_Area[1][bin], JetPt_R_Area[1][bin], JetPt_Area_Wide[1][bin], Z_Area_Wide[1][bin], ZPt_Area_Wide[1][bin], JetPt_R_Area_Wide[1][bin], JetPtAreaVsJetPtCS[1][bin], NULL, JetAreaVsJetNConstArea[1][bin], JetAreaVsJetNConstCS[1][bin], JetPtAreaVsJetNConstArea[1][bin], JetPtCSVsJetNConstCS[1][bin]);
    //     sumw2 += fillHistos(Form("./sWeights/NewsWeight_PT_%i_%i_Cen_10_20.root", (int)ptlow, (int)pthigh), mode, D0PtUnweighted[2][bin], D0Pt[2][bin], JetPt[2][bin], Z[2][bin], R[2][bin], ZPt[2][bin], 10, 20, hArea_2, JetEta_2, D0Eta_2, JetPhi_2, D0Phi_2, D0JetEta_2, Eta_2, JetPt_R[2][bin], JetPt_Z[2][bin], widebin, JetPt_Area[2][bin], Z_Area[2][bin], ZPt_Area[2][bin], JetPt_R_Area[2][bin], JetPt_Area_Wide[2][bin], Z_Area_Wide[2][bin], ZPt_Area_Wide[2][bin], JetPt_R_Area_Wide[2][bin], JetPtAreaVsJetPtCS[2][bin], NULL, JetAreaVsJetNConstArea[2][bin], JetAreaVsJetNConstCS[2][bin], JetPtAreaVsJetNConstArea[2][bin], JetPtCSVsJetNConstCS[2][bin]);
    //     sumw3 += fillHistos(Form("./sWeights/NewsWeight_PT_%i_%i_Cen_20_30.root", (int)ptlow, (int)pthigh), mode, D0PtUnweighted[3][bin], D0Pt[3][bin], JetPt[3][bin], Z[3][bin], R[3][bin], ZPt[3][bin], 20, 30, hArea_3, JetEta_3, D0Eta_3, JetPhi_3, D0Phi_3, D0JetEta_3, Eta_3, JetPt_R[3][bin], JetPt_Z[3][bin], widebin, JetPt_Area[3][bin], Z_Area[3][bin], ZPt_Area[3][bin], JetPt_R_Area[3][bin], JetPt_Area_Wide[3][bin], Z_Area_Wide[3][bin], ZPt_Area_Wide[3][bin], JetPt_R_Area_Wide[3][bin], JetPtAreaVsJetPtCS[3][bin], NULL, JetAreaVsJetNConstArea[3][bin], JetAreaVsJetNConstCS[3][bin], JetPtAreaVsJetNConstArea[3][bin], JetPtCSVsJetNConstCS[3][bin]);
    //     sumw4 += fillHistos(Form("./sWeights/NewsWeight_PT_%i_%i_Cen_30_40.root", (int)ptlow, (int)pthigh), mode, D0PtUnweighted[4][bin], D0Pt[4][bin], JetPt[4][bin], Z[4][bin], R[4][bin], ZPt[4][bin], 30, 40, hArea_4, JetEta_4, D0Eta_4, JetPhi_4, D0Phi_4, D0JetEta_4, Eta_4, JetPt_R[4][bin], JetPt_Z[4][bin], widebin, JetPt_Area[4][bin], Z_Area[4][bin], ZPt_Area[4][bin], JetPt_R_Area[4][bin], JetPt_Area_Wide[4][bin], Z_Area_Wide[4][bin], ZPt_Area_Wide[4][bin], JetPt_R_Area_Wide[4][bin], JetPtAreaVsJetPtCS[4][bin], NULL, JetAreaVsJetNConstArea[4][bin], JetAreaVsJetNConstCS[4][bin], JetPtAreaVsJetNConstArea[4][bin], JetPtCSVsJetNConstCS[4][bin]);
    //     sumw6 += fillHistos(Form("./sWeights/NewsWeight_PT_%i_%i_Cen_10_20.root", (int)ptlow, (int)pthigh), mode, D0PtUnweighted[6][bin], D0Pt[6][bin], JetPt[6][bin], Z[6][bin], R[6][bin], ZPt[6][bin], 10, 20, hArea_6, JetEta_6, D0Eta_6, JetPhi_6, D0Phi_6, D0JetEta_6, Eta_6, JetPt_R[6][bin], JetPt_Z[6][bin], widebin, JetPt_Area[6][bin], Z_Area[6][bin], ZPt_Area[6][bin], JetPt_R_Area[6][bin], JetPt_Area_Wide[6][bin], Z_Area_Wide[6][bin], ZPt_Area_Wide[6][bin], JetPt_R_Area_Wide[6][bin], JetPtAreaVsJetPtCS[6][bin], NULL, JetAreaVsJetNConstArea[6][bin], JetAreaVsJetNConstCS[6][bin], JetPtAreaVsJetNConstArea[6][bin], JetPtCSVsJetNConstCS[6][bin]);
    //     sumw6 += fillHistos(Form("./sWeights/NewsWeight_PT_%i_%i_Cen_20_30.root", (int)ptlow, (int)pthigh), mode, D0PtUnweighted[6][bin], D0Pt[6][bin], JetPt[6][bin], Z[6][bin], R[6][bin], ZPt[6][bin], 20, 30, hArea_6, JetEta_6, D0Eta_6, JetPhi_6, D0Phi_6, D0JetEta_6, Eta_6, JetPt_R[6][bin], JetPt_Z[6][bin], widebin, JetPt_Area[6][bin], Z_Area[6][bin], ZPt_Area[6][bin], JetPt_R_Area[6][bin], JetPt_Area_Wide[6][bin], Z_Area_Wide[6][bin], ZPt_Area_Wide[6][bin], JetPt_R_Area_Wide[6][bin], JetPtAreaVsJetPtCS[6][bin], NULL, JetAreaVsJetNConstArea[6][bin], JetAreaVsJetNConstCS[6][bin], JetPtAreaVsJetNConstArea[6][bin], JetPtCSVsJetNConstCS[6][bin]);
    //     sumw6 += fillHistos(Form("./sWeights/NewsWeight_PT_%i_%i_Cen_30_40.root", (int)ptlow, (int)pthigh), mode, D0PtUnweighted[6][bin], D0Pt[6][bin], JetPt[6][bin], Z[6][bin], R[6][bin], ZPt[6][bin], 30, 40, hArea_6, JetEta_6, D0Eta_6, JetPhi_6, D0Phi_6, D0JetEta_6, Eta_6, JetPt_R[6][bin], JetPt_Z[6][bin], widebin, JetPt_Area[6][bin], Z_Area[6][bin], ZPt_Area[6][bin], JetPt_R_Area[6][bin], JetPt_Area_Wide[6][bin], Z_Area_Wide[6][bin], ZPt_Area_Wide[6][bin], JetPt_R_Area_Wide[6][bin], JetPtAreaVsJetPtCS[6][bin], NULL, JetAreaVsJetNConstArea[6][bin], JetAreaVsJetNConstCS[6][bin], JetPtAreaVsJetNConstArea[6][bin], JetPtCSVsJetNConstCS[6][bin]);
    //     sumw7 += fillHistos(Form("./sWeights/NewsWeight_PT_%i_%i_Cen_40_80.root", (int)ptlow, (int)pthigh), mode, D0PtUnweighted[7][bin], D0Pt[7][bin], JetPt[7][bin], Z[7][bin], R[7][bin], ZPt[7][bin], 40, 80, hArea_7, JetEta_7, D0Eta_7, JetPhi_7, D0Phi_7, D0JetEta_7, Eta_7, JetPt_R[7][bin], JetPt_Z[7][bin], widebin, JetPt_Area[7][bin], Z_Area[7][bin], ZPt_Area[7][bin], JetPt_R_Area[7][bin], JetPt_Area_Wide[7][bin], Z_Area_Wide[7][bin], ZPt_Area_Wide[7][bin], JetPt_R_Area_Wide[7][bin], JetPtAreaVsJetPtCS[7][bin], NULL, JetAreaVsJetNConstArea[7][bin], JetAreaVsJetNConstCS[7][bin], JetPtAreaVsJetNConstArea[7][bin], JetPtCSVsJetNConstCS[7][bin]);

    //     cout << "Sum = " << sumw1 << "\t" << sumw6 << "\t" << sumw7 << endl;
    //     cout << "Sum From hist = " << D0Pt[1][0]->Integral() << "\t" << D0Pt[6][0]->Integral() << "\t" << D0Pt[7][0]->Integral() << endl;
    
    // }

    double tmpval = fillHistos(Form("./sWeights/NewsWeight_PT_%i_%i_Cen_0_80.root",  D0PtBin[0], D0PtBin[nd0ptbins]), mode, SYSMODE, NULL, NULL, NULL, NULL, NULL, NULL, 0, 80, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, widebin, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, MuTreeForMultiFold);

    for (int bin = 0; bin < nd0ptbins; bin++){
        sumw  += fillHistos(Form("./sWeights/NewsWeight_PT_%i_%i_Cen_0_80.root",  D0PtBin[bin], D0PtBin[bin+1]), mode, SYSMODE, D0PtUnweighted[0][bin], D0Pt[0][bin], JetPt[0][bin], Z[0][bin], R[0][bin], ZPt[0][bin], 0, 80, hArea, JetEta, D0Eta, JetPhi, D0Phi, D0JetEta, Eta, JetPt_R[0][bin], JetPt_Z[0][bin], widebin, JetPt_Area[0][bin], Z_Area[0][bin], ZPt_Area[0][bin], JetPt_R_Area[0][bin], JetPt_Area_Wide[0][bin], Z_Area_Wide[0][bin], ZPt_Area_Wide[0][bin], JetPt_R_Area_Wide[0][bin], JetPtAreaVsJetPtCS[0][bin], centrality, JetAreaVsJetNConstArea[0][bin], JetAreaVsJetNConstCS[0][bin], JetPtAreaVsJetNConstArea[0][bin], JetPtCSVsJetNConstCS[0][bin], JetCentVsJetNConst[0][bin], JetRefMult[0][bin]);
        sumw1 += fillHistos(Form("./sWeights/NewsWeight_PT_%i_%i_Cen_0_10.root",  D0PtBin[bin], D0PtBin[bin+1]), mode, SYSMODE, D0PtUnweighted[1][bin], D0Pt[1][bin], JetPt[1][bin], Z[1][bin], R[1][bin], ZPt[1][bin], 0, 10, hArea_1, JetEta_1, D0Eta_1, JetPhi_1, D0Phi_1, D0JetEta_1, Eta_1, JetPt_R[1][bin], JetPt_Z[1][bin], widebin, JetPt_Area[1][bin], Z_Area[1][bin], ZPt_Area[1][bin], JetPt_R_Area[1][bin], JetPt_Area_Wide[1][bin], Z_Area_Wide[1][bin], ZPt_Area_Wide[1][bin], JetPt_R_Area_Wide[1][bin], JetPtAreaVsJetPtCS[1][bin], NULL, JetAreaVsJetNConstArea[1][bin], JetAreaVsJetNConstCS[1][bin], JetPtAreaVsJetNConstArea[1][bin], JetPtCSVsJetNConstCS[1][bin], NULL, JetRefMult[1][bin]);
        sumw2 += fillHistos(Form("./sWeights/NewsWeight_PT_%i_%i_Cen_10_20.root", D0PtBin[bin], D0PtBin[bin+1]), mode, SYSMODE, D0PtUnweighted[2][bin], D0Pt[2][bin], JetPt[2][bin], Z[2][bin], R[2][bin], ZPt[2][bin], 10, 20, hArea_2, JetEta_2, D0Eta_2, JetPhi_2, D0Phi_2, D0JetEta_2, Eta_2, JetPt_R[2][bin], JetPt_Z[2][bin], widebin, JetPt_Area[2][bin], Z_Area[2][bin], ZPt_Area[2][bin], JetPt_R_Area[2][bin], JetPt_Area_Wide[2][bin], Z_Area_Wide[2][bin], ZPt_Area_Wide[2][bin], JetPt_R_Area_Wide[2][bin], JetPtAreaVsJetPtCS[2][bin], NULL, JetAreaVsJetNConstArea[2][bin], JetAreaVsJetNConstCS[2][bin], JetPtAreaVsJetNConstArea[2][bin], JetPtCSVsJetNConstCS[2][bin], NULL, JetRefMult[2][bin]);
        sumw3 += fillHistos(Form("./sWeights/NewsWeight_PT_%i_%i_Cen_20_30.root", D0PtBin[bin], D0PtBin[bin+1]), mode, SYSMODE, D0PtUnweighted[3][bin], D0Pt[3][bin], JetPt[3][bin], Z[3][bin], R[3][bin], ZPt[3][bin], 20, 30, hArea_3, JetEta_3, D0Eta_3, JetPhi_3, D0Phi_3, D0JetEta_3, Eta_3, JetPt_R[3][bin], JetPt_Z[3][bin], widebin, JetPt_Area[3][bin], Z_Area[3][bin], ZPt_Area[3][bin], JetPt_R_Area[3][bin], JetPt_Area_Wide[3][bin], Z_Area_Wide[3][bin], ZPt_Area_Wide[3][bin], JetPt_R_Area_Wide[3][bin], JetPtAreaVsJetPtCS[3][bin], NULL, JetAreaVsJetNConstArea[3][bin], JetAreaVsJetNConstCS[3][bin], JetPtAreaVsJetNConstArea[3][bin], JetPtCSVsJetNConstCS[3][bin], NULL, JetRefMult[3][bin]);
        sumw4 += fillHistos(Form("./sWeights/NewsWeight_PT_%i_%i_Cen_30_40.root", D0PtBin[bin], D0PtBin[bin+1]), mode, SYSMODE, D0PtUnweighted[4][bin], D0Pt[4][bin], JetPt[4][bin], Z[4][bin], R[4][bin], ZPt[4][bin], 30, 40, hArea_4, JetEta_4, D0Eta_4, JetPhi_4, D0Phi_4, D0JetEta_4, Eta_4, JetPt_R[4][bin], JetPt_Z[4][bin], widebin, JetPt_Area[4][bin], Z_Area[4][bin], ZPt_Area[4][bin], JetPt_R_Area[4][bin], JetPt_Area_Wide[4][bin], Z_Area_Wide[4][bin], ZPt_Area_Wide[4][bin], JetPt_R_Area_Wide[4][bin], JetPtAreaVsJetPtCS[4][bin], NULL, JetAreaVsJetNConstArea[4][bin], JetAreaVsJetNConstCS[4][bin], JetPtAreaVsJetNConstArea[4][bin], JetPtCSVsJetNConstCS[4][bin], NULL, JetRefMult[4][bin]);
        sumw6 += fillHistos(Form("./sWeights/NewsWeight_PT_%i_%i_Cen_10_20.root", D0PtBin[bin], D0PtBin[bin+1]), mode, SYSMODE, D0PtUnweighted[6][bin], D0Pt[6][bin], JetPt[6][bin], Z[6][bin], R[6][bin], ZPt[6][bin], 10, 20, hArea_6, JetEta_6, D0Eta_6, JetPhi_6, D0Phi_6, D0JetEta_6, Eta_6, JetPt_R[6][bin], JetPt_Z[6][bin], widebin, JetPt_Area[6][bin], Z_Area[6][bin], ZPt_Area[6][bin], JetPt_R_Area[6][bin], JetPt_Area_Wide[6][bin], Z_Area_Wide[6][bin], ZPt_Area_Wide[6][bin], JetPt_R_Area_Wide[6][bin], JetPtAreaVsJetPtCS[6][bin], NULL, JetAreaVsJetNConstArea[6][bin], JetAreaVsJetNConstCS[6][bin], JetPtAreaVsJetNConstArea[6][bin], JetPtCSVsJetNConstCS[6][bin], NULL, JetRefMult[6][bin]);
        sumw6 += fillHistos(Form("./sWeights/NewsWeight_PT_%i_%i_Cen_20_30.root", D0PtBin[bin], D0PtBin[bin+1]), mode, SYSMODE, D0PtUnweighted[6][bin], D0Pt[6][bin], JetPt[6][bin], Z[6][bin], R[6][bin], ZPt[6][bin], 20, 30, hArea_6, JetEta_6, D0Eta_6, JetPhi_6, D0Phi_6, D0JetEta_6, Eta_6, JetPt_R[6][bin], JetPt_Z[6][bin], widebin, JetPt_Area[6][bin], Z_Area[6][bin], ZPt_Area[6][bin], JetPt_R_Area[6][bin], JetPt_Area_Wide[6][bin], Z_Area_Wide[6][bin], ZPt_Area_Wide[6][bin], JetPt_R_Area_Wide[6][bin], JetPtAreaVsJetPtCS[6][bin], NULL, JetAreaVsJetNConstArea[6][bin], JetAreaVsJetNConstCS[6][bin], JetPtAreaVsJetNConstArea[6][bin], JetPtCSVsJetNConstCS[6][bin], NULL, JetRefMult[6][bin]);
        sumw6 += fillHistos(Form("./sWeights/NewsWeight_PT_%i_%i_Cen_30_40.root", D0PtBin[bin], D0PtBin[bin+1]), mode, SYSMODE, D0PtUnweighted[6][bin], D0Pt[6][bin], JetPt[6][bin], Z[6][bin], R[6][bin], ZPt[6][bin], 30, 40, hArea_6, JetEta_6, D0Eta_6, JetPhi_6, D0Phi_6, D0JetEta_6, Eta_6, JetPt_R[6][bin], JetPt_Z[6][bin], widebin, JetPt_Area[6][bin], Z_Area[6][bin], ZPt_Area[6][bin], JetPt_R_Area[6][bin], JetPt_Area_Wide[6][bin], Z_Area_Wide[6][bin], ZPt_Area_Wide[6][bin], JetPt_R_Area_Wide[6][bin], JetPtAreaVsJetPtCS[6][bin], NULL, JetAreaVsJetNConstArea[6][bin], JetAreaVsJetNConstCS[6][bin], JetPtAreaVsJetNConstArea[6][bin], JetPtCSVsJetNConstCS[6][bin], NULL, JetRefMult[6][bin]);
        sumw7 += fillHistos(Form("./sWeights/NewsWeight_PT_%i_%i_Cen_40_80.root", D0PtBin[bin], D0PtBin[bin+1]), mode, SYSMODE, D0PtUnweighted[7][bin], D0Pt[7][bin], JetPt[7][bin], Z[7][bin], R[7][bin], ZPt[7][bin], 40, 80, hArea_7, JetEta_7, D0Eta_7, JetPhi_7, D0Phi_7, D0JetEta_7, Eta_7, JetPt_R[7][bin], JetPt_Z[7][bin], widebin, JetPt_Area[7][bin], Z_Area[7][bin], ZPt_Area[7][bin], JetPt_R_Area[7][bin], JetPt_Area_Wide[7][bin], Z_Area_Wide[7][bin], ZPt_Area_Wide[7][bin], JetPt_R_Area_Wide[7][bin], JetPtAreaVsJetPtCS[7][bin], NULL, JetAreaVsJetNConstArea[7][bin], JetAreaVsJetNConstCS[7][bin], JetPtAreaVsJetNConstArea[7][bin], JetPtCSVsJetNConstCS[7][bin], NULL, JetRefMult[7][bin]);

        cout << "Sum = " << sumw1 << "\t" << sumw6 << "\t" << sumw7 << endl;
        cout << "Sum From hist = " << JetPt_R_Area_Wide[1][0]->Integral() << "\t" << JetPt_R_Area_Wide[6][0]->Integral() << "\t" << JetPt_R_Area_Wide[7][0]->Integral() << endl;
    }

    // Remove negative yield bins for the array of histograms I declared earlier.

    
    for (int i = 0; i < 8; i++)
    { // Centrality Bins
        for (int bin = 0; bin < nd0ptbins; bin++)
        { // D0 pT Bins
            RemoveNegativeYieldBins(D0PtUnweighted[i][bin]);
            RemoveNegativeYieldBins(D0Pt[i][bin]);
            RemoveNegativeYieldBins(JetPt[i][bin]);
            RemoveNegativeYieldBins(Z[i][bin]);
            RemoveNegativeYieldBins(R[i][bin]);
            RemoveNegativeYieldBins(ZPt[i][bin]);
            RemoveNegativeYieldBins(ZPt_Area[i][bin]);
            RemoveNegativeYieldBins(ZPt_Area_Wide[i][bin]);
            RemoveNegativeYieldBins(JetPt_Area_Wide[i][bin]);
            RemoveNegativeYieldBins(JetPt_R[i][bin]);
            // RemoveNegativeYieldBins(JetPt_Z[i][bin]);
            RemoveNegativeYieldBins(JetCentVsJetNConst[i][bin]);
            RemoveNegativeYieldBins(JetRefMult[i][bin]);
            RemoveNegativeYieldBins(JetPt_R_Area_Wide[i][bin]);
        }
    }

    // Now I want to add all the D0 pT bins to the first D0 pT bin for each centrality bin.

    for (int i = 0; i < 8; i++){ // Centrality Bins
        for (int bin = 1; bin < nd0ptbins; bin++)
        { // D0 pT Bins
            D0PtUnweighted[i][0]->Add(D0PtUnweighted[i][bin]);
            D0Pt[i][0]->Add(D0Pt[i][bin]);
            JetPt[i][0]->Add(JetPt[i][bin]);
            Z[i][0]->Add(Z[i][bin]);
            R[i][0]->Add(R[i][bin]);
            ZPt[i][0]->Add(ZPt[i][bin]);
            ZPt_Area[i][0]->Add(ZPt_Area[i][bin]);
            ZPt_Area_Wide[i][0]->Add(ZPt_Area_Wide[i][bin]);
            JetPt_Area_Wide[i][0]->Add(JetPt_Area_Wide[i][bin]);
            JetPt_R[i][0]->Add(JetPt_R[i][bin]);
            JetPt_Z[i][0]->Add(JetPt_Z[i][bin]);
            JetPt_R_Area_Wide[i][0]->Add(JetPt_R_Area_Wide[i][bin]);
            JetRefMult[i][0]->Add(JetRefMult[i][bin]);
        }
    }

    if (!save)
    {
        DN(D0Pt[0][0], 14, 4, sumw);
        DN(D0Pt[1][0], 2, 21, sumw1);
        DN(D0Pt[2][0], 9, 22, sumw2);
        DN(D0Pt[3][0], 32, 23, sumw3);
        DN(D0Pt[4][0], 6, 24, sumw4);
        DN(D0Pt[5][0], 1, 8, sumw5);

        DN(JetPt[0][0], 14, 4, 1);  // sumw);
        DN(JetPt[1][0], 2, 21, 1);  // sumw1);
        DN(JetPt[2][0], 9, 22, 1);  // sumw2);
        DN(JetPt[3][0], 32, 23, 1); // sumw3);
        DN(JetPt[4][0], 6, 24, 1);  // sumw4);
        DN(JetPt[5][0], 1, 8, 1);   // sumw5);

        DN(R[0][0], 14, 4, sumw);
        DN(R[1][0], 2, 21, sumw1);
        DN(R[2][0], 9, 22, sumw2);
        DN(R[3][0], 32, 23, sumw3);
        DN(R[4][0], 6, 24, sumw4);
        DN(R[5][0], 1, 8, sumw5);

        DN(Z[0][0], 14, 4, sumw);
        DN(Z[1][0], 2, 21, sumw1);
        DN(Z[2][0], 9, 22, sumw2);
        DN(Z[3][0], 32, 23, sumw3);
        DN(Z[4][0], 6, 24, sumw4);
        DN(Z[5][0], 1, 8, sumw5);
    }

    D0Pt[0][0]->GetXaxis()->SetTitle("D^{0} p_{T} (GeV/#it{c})");
    JetPt[0][0]->GetXaxis()->SetTitle("Jet p_{T} (GeV/#it{c})");
    Z[0][0]->GetXaxis()->SetTitle("#it{z}");
    R[0][0]->GetXaxis()->SetTitle("#DeltaR");

    D0Pt[0][0]->GetYaxis()->SetTitle("dN/d#it{p}_{T} (#it{c}/GeV)");
    JetPt[0][0]->GetYaxis()->SetTitle("dN/d#it{p}_{T} (#it{c}/GeV)");
    Z[0][0]->GetYaxis()->SetTitle("(1/N_{jet}) dN/d#it{z}");
    R[0][0]->GetYaxis()->SetTitle("(1/N_{jet}) ddN/d#DeltaR");

    D0Pt[0][0]->GetXaxis()->SetRangeUser(0, 10);
    JetPt[0][0]->GetXaxis()->SetRangeUser(3, 40);
    R[0][0]->GetXaxis()->SetRangeUser(0, 0.4);
    Z[0][0]->GetXaxis()->SetRangeUser(0, 1);

    D0Pt[0][0]->GetYaxis()->SetRangeUser(0.00001, D0Pt[0][0]->GetMaximum() * 10.5);
    JetPt[0][0]->GetYaxis()->SetRangeUser(0.00001, JetPt[0][0]->GetMaximum() * 10.5);
    R[0][0]->GetYaxis()->SetRangeUser(0.1, R[0][0]->GetMaximum() * 10.5);
    Z[0][0]->GetYaxis()->SetRangeUser(0.001, Z[0][0]->GetMaximum() * 10.5);

    if (!SYS)
    {
        TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
        leg->SetTextSize(0.04);
        leg->AddEntry(R[0][0], "0-80%", "PE");
        leg->AddEntry(R[1][0], "0-10%", "PE");
        leg->AddEntry(R[2][0], "10-20%", "PE");
        leg->AddEntry(R[3][0], "20-30%", "PE");
        leg->AddEntry(R[4][0], "30-40%", "PE");
        leg->AddEntry(R[5][0], "40-80%", "PE");
        TCanvas *c11 = new TCanvas("c11", "c11");
        D0Pt[0][0]->DrawClone("PE");
        D0Pt[1][0]->DrawClone("PE same");
        D0Pt[2][0]->DrawClone("PE same");
        D0Pt[3][0]->DrawClone("PE same");
        D0Pt[4][0]->DrawClone("PE same");
        D0Pt[5][0]->DrawClone("PE same");
        leg->Draw("same");
        gPad->SetLogy();

        TCanvas *c21 = new TCanvas("c21", "c21");
        JetPt[0][0]->DrawClone("PE");
        // TCanvas *c211 =new TCanvas("c211","c211");
        JetPt[1][0]->DrawClone("PE same");
        // TCanvas *c212 =new TCanvas("c212","c212");
        JetPt[2][0]->DrawClone("PE same");
        // TCanvas *c213 =new TCanvas("c213","c213");
        JetPt[3][0]->DrawClone("PE same");
        // TCanvas *c21 =new TCanvas("c24","c24");
        JetPt[4][0]->DrawClone("PE same");
        JetPt[5][0]->DrawClone("PE same");
        leg->Draw("same");
        gPad->SetLogy();
        TCanvas *c31 = new TCanvas("c31", "c31");
        Z[0][0]->DrawClone("PE");
        Z[1][0]->DrawClone("PE same");
        Z[2][0]->DrawClone("PE same");
        Z[3][0]->DrawClone("PE same");
        Z[4][0]->DrawClone("PE same");
        Z[5][0]->DrawClone("PE same");
        leg->Draw("same");
        gPad->SetLogy();
        TCanvas *c41 = new TCanvas("c41", "c41");
        R[0][0]->DrawClone("PE");
        R[1][0]->DrawClone("PE same");
        R[2][0]->DrawClone("PE same");
        R[3][0]->DrawClone("PE same");
        R[4][0]->DrawClone("PE same");
        R[5][0]->DrawClone("PE same");
        leg->Draw("same");
        gPad->SetLogy();

        PT2D->GetXaxis()->SetTitle("Jet p^{corr.}_{T} (GeV/c)");
        PT2D->GetYaxis()->SetTitle("D^{0} p_{T} (GeV/c)");
        PT2D->GetZaxis()->SetRangeUser(2, PT2D->GetMaximum() * 5);
        TCanvas *c411 = new TCanvas("c411", "c411");
        PT2D->Draw("COLZ");
    }


    if (save)
    {
        outformultifold->cd();
        MuTreeForMultiFold->Write();
        outformultifold->Close();

        TFile *outf = new TFile(NAME.Data(), "RECREATE");

        cout << "File Name = " << NAME.Data() << endl;
        centrality->Write();

        cout << "Mode ==== " << mode << endl;
        cout << "Integral here " << D0Pt[1][0]->Integral() << " " << sumw1 << endl;
        cout << "Integral here " << D0Pt[6][0]->Integral() << " " << sumw6 << endl;
        cout << "Integral here " << D0Pt[7][0]->Integral() << " " << sumw7 << endl;
        cout << "===================================" << endl;
        cout << "Integral here " << ZPt[1][0]->Integral() << " " << sumw1 << endl;
        cout << "Integral here " << ZPt[6][0]->Integral() << " " << sumw6 << endl;
        cout << "Integral here " << ZPt[7][0]->Integral() << " " << sumw7 << endl;
        cout << "===================================" << endl;
        cout << "Integral here " << ZPt_Area_Wide[1][0]->Integral(0, ZPt_Area_Wide[1][0]->GetNbinsX()+1, 0, ZPt_Area_Wide[1][0]->GetNbinsY()+1) << " " << sumw1 << endl;
        cout << "Integral here " << ZPt_Area_Wide[6][0]->Integral(0, ZPt_Area_Wide[6][0]->GetNbinsX()+1, 0, ZPt_Area_Wide[6][0]->GetNbinsY()+1) << " " << sumw6 << endl;
        cout << "Integral here " << ZPt_Area_Wide[7][0]->Integral(0, ZPt_Area_Wide[7][0]->GetNbinsX()+1, 0, ZPt_Area_Wide[7][0]->GetNbinsY()+1) << " " << sumw7 << endl;

        D0PtUnweighted[0][0]->Write("D0PtUnweighted");
        D0PtUnweighted[1][0]->Write("D0PtUnweighted_0_10");
        D0PtUnweighted[2][0]->Write("D0PtUnweighted_10_20");
        D0PtUnweighted[3][0]->Write("D0PtUnweighted_20_30");
        D0PtUnweighted[4][0]->Write("D0PtUnweighted_30_40");
        D0PtUnweighted[6][0]->Write("D0PtUnweighted_10_40");
        D0PtUnweighted[7][0]->Write("D0PtUnweighted_40_80");

        D0Pt[0][0]->Write("D0Pt");
        D0Pt[1][0]->Write("D0Pt_0_10");
        D0Pt[2][0]->Write("D0Pt_10_20");
        D0Pt[3][0]->Write("D0Pt_20_30");
        D0Pt[4][0]->Write("D0Pt_30_40");
        D0Pt[6][0]->Write("D0Pt_10_40");
        D0Pt[7][0]->Write("D0Pt_40_80");

        D0JPt_5->Write("D0JPt_40_80");

        JetPt[0][0]->Write("JetPt");
        JetPt[1][0]->Write("JetPt_0_10");
        JetPt[2][0]->Write("JetPt_10_20");
        JetPt[3][0]->Write("JetPt_20_30");
        JetPt[4][0]->Write("JetPt_30_40");
        JetPt[6][0]->Write("JetPt_10_40");
        JetPt[7][0]->Write("JetPt_40_80");

        Z[0][0]->Write("Z");
        Z[1][0]->Write("Z_0_10");
        Z[2][0]->Write("Z_10_20");
        Z[3][0]->Write("Z_20_30");
        Z[4][0]->Write("Z_30_40");
        Z[6][0]->Write("Z_10_40");
        Z[7][0]->Write("Z_40_80");

        // cout << "Negative Yield Bins " << ZPt_1->GetBinContent(ZPt_1->GetMinimumBin()) << endl;
        // cout << "Negative Yield Bins " << ZPt_6->GetBinContent(ZPt_6->GetMinimumBin()) << endl;
        // cout << "Negative Yield Bins " << ZPt_7->GetBinContent(ZPt_7->GetMinimumBin()) << endl;

        ZPt[0][0]->Write("ZPt");
        ZPt[1][0]->Write("ZPt_0_10");
        ZPt[2][0]->Write("ZPt_10_20");
        ZPt[3][0]->Write("ZPt_20_30");
        ZPt[4][0]->Write("ZPt_30_40");
        ZPt[6][0]->Write("ZPt_10_40");
        ZPt[7][0]->Write("ZPt_40_80");

        ZPt[0][0]->ProjectionX()->Write("ZPt_Px");
        ZPt[1][0]->ProjectionX()->Write("ZPt_Px_0_10");
        ZPt[2][0]->ProjectionX()->Write("ZPt_Px_10_20");
        ZPt[3][0]->ProjectionX()->Write("ZPt_Px_20_30");
        ZPt[4][0]->ProjectionX()->Write("ZPt_Px_30_40");
        ZPt[6][0]->ProjectionX()->Write("ZPt_Px_10_40");
        ZPt[7][0]->ProjectionX()->Write("ZPt_Px_40_80");

        ZPt[0][0]->ProjectionY()->Write("ZPt_Py");
        ZPt[1][0]->ProjectionY()->Write("ZPt_Py_0_10");
        ZPt[2][0]->ProjectionY()->Write("ZPt_Py_10_20");
        ZPt[3][0]->ProjectionY()->Write("ZPt_Py_20_30");
        ZPt[4][0]->ProjectionY()->Write("ZPt_Py_30_40");
        ZPt[6][0]->ProjectionY()->Write("ZPt_Py_10_40");
        ZPt[7][0]->ProjectionY()->Write("ZPt_Py_40_80");

        R[0][0]->Write("R");
        R[1][0]->Write("R_0_10");
        R[2][0]->Write("R_10_20");
        R[3][0]->Write("R_20_30");
        R[4][0]->Write("R_30_40");
        R[6][0]->Write("R_10_40");
        R[7][0]->Write("R_40_80");

        hArea->Write("JetArea_0_80");
        Hard2D2_5->Write("Hard2D2_40_80");

        JetEta->Write("JetEta");
        JetEta_1->Write("JetEta_1");
        JetEta_2->Write("JetEta_2");
        JetEta_3->Write("JetEta_3");
        JetEta_4->Write("JetEta_4");
        JetEta_5->Write("JetEta_5");
        JetEta_6->Write("JetEta_6");
        JetEta_7->Write("JetEta_7");
        D0Eta->Write("D0Eta");
        D0Eta_1->Write("D0Eta_1");
        D0Eta_2->Write("D0Eta_2");
        D0Eta_3->Write("D0Eta_3");
        D0Eta_4->Write("D0Eta_4");
        D0Eta_5->Write("D0Eta_5");
        D0Eta_6->Write("D0Eta_6");
        D0Eta_7->Write("D0Eta_7");

        Eta->Write("Eta");
        Eta_1->Write("Eta_1");
        Eta_2->Write("Eta_2");
        Eta_3->Write("Eta_3");
        Eta_4->Write("Eta_4");
        Eta_5->Write("Eta_5");
        Eta_6->Write("Eta_6");
        Eta_7->Write("Eta_7");

        D0JetEta->Write("D0JetEta");
        D0JetEta_1->Write("D0JetEta_1");
        D0JetEta_2->Write("D0JetEta_2");
        D0JetEta_3->Write("D0JetEta_3");
        D0JetEta_4->Write("D0JetEta_4");
        D0JetEta_5->Write("D0JetEta_5");
        D0JetEta_6->Write("D0JetEta_6");
        D0JetEta_7->Write("D0JetEta_7");

        hArea_1->Write("hArea_1");
        hArea_2->Write("hArea_2");
        hArea_3->Write("hArea_3");
        hArea_4->Write("hArea_4");
        hArea_5->Write("hArea_5");
        hArea_6->Write("hArea_6");
        hArea_7->Write("hArea_7");
        JetPhi->Write("JetPhi");
        JetPhi_1->Write("JetPhi_1");
        JetPhi_2->Write("JetPhi_2");
        JetPhi_3->Write("JetPhi_3");
        JetPhi_4->Write("JetPhi_4");
        JetPhi_5->Write("JetPhi_5");
        JetPhi_6->Write("JetPhi_6");
        JetPhi_7->Write("JetPhi_7");

        D0Phi->Write("D0Phi");
        D0Phi_1->Write("D0Phi_1");
        D0Phi_2->Write("D0Phi_2");
        D0Phi_3->Write("D0Phi_3");
        D0Phi_4->Write("D0Phi_4");
        D0Phi_5->Write("D0Phi_5");
        D0Phi_6->Write("D0Phi_6");
        D0Phi_7->Write("D0Phi_7");

        JetPt_R[0][0]->Write("JetPt_R");
        JetPt_R[1][0]->Write("JetPt_R_0_10");
        JetPt_R[2][0]->Write("JetPt_R_10_20");
        JetPt_R[3][0]->Write("JetPt_R_20_30");
        JetPt_R[4][0]->Write("JetPt_R_30_40");
        JetPt_R[6][0]->Write("JetPt_R_10_40");
        JetPt_R[7][0]->Write("JetPt_R_40_80");

        JetPt_Z[0][0]->Write("JetPt_Z_0_80");
        JetPt_Z[1][0]->Write("JetPt_Z_0_10");
        JetPt_Z[6][0]->Write("JetPt_Z_10_40");
        JetPt_Z[7][0]->Write("JetPt_Z_40_80");
        
        if (mode == "CS"){
            JetPt_Area[0][0]->Write("JetPt_Area");
            JetPt_Area[1][0]->Write("JetPt_Area_0_10");
            JetPt_Area[2][0]->Write("JetPt_Area_10_20");
            JetPt_Area[3][0]->Write("JetPt_Area_20_30");
            JetPt_Area[4][0]->Write("JetPt_Area_30_40");
            JetPt_Area[6][0]->Write("JetPt_Area_10_40");
            JetPt_Area[7][0]->Write("JetPt_Area_40_80");

            Z_Area[0][0]->Write("Z_Area");
            Z_Area[1][0]->Write("Z_Area_0_10");
            Z_Area[2][0]->Write("Z_Area_10_20");
            Z_Area[3][0]->Write("Z_Area_20_30");
            Z_Area[4][0]->Write("Z_Area_30_40");
            Z_Area[6][0]->Write("Z_Area_10_40");
            Z_Area[7][0]->Write("Z_Area_40_80");

            ZPt_Area[0][0]->Write("ZPt_Area");
            ZPt_Area[1][0]->Write("ZPt_Area_0_10");
            ZPt_Area[2][0]->Write("ZPt_Area_10_20");
            ZPt_Area[3][0]->Write("ZPt_Area_20_30");
            ZPt_Area[4][0]->Write("ZPt_Area_30_40");
            ZPt_Area[6][0]->Write("ZPt_Area_10_40");
            ZPt_Area[7][0]->Write("ZPt_Area_40_80");

            JetPt_R_Area[0][0]->Write("JetPt_R_Area");
            JetPt_R_Area[1][0]->Write("JetPt_R_Area_0_10");
            JetPt_R_Area[2][0]->Write("JetPt_R_Area_10_20");
            JetPt_R_Area[3][0]->Write("JetPt_R_Area_20_30");
            JetPt_R_Area[4][0]->Write("JetPt_R_Area_30_40");
            JetPt_R_Area[6][0]->Write("JetPt_R_Area_10_40");
            JetPt_R_Area[7][0]->Write("JetPt_R_Area_40_80");

            JetPt_Area_Wide[0][0]->Write("JetPt_Area_Wide");
            JetPt_Area_Wide[1][0]->Write("JetPt_Area_Wide_0_10");
            JetPt_Area_Wide[2][0]->Write("JetPt_Area_Wide_10_20");
            JetPt_Area_Wide[3][0]->Write("JetPt_Area_Wide_20_30");
            JetPt_Area_Wide[4][0]->Write("JetPt_Area_Wide_30_40");
            JetPt_Area_Wide[6][0]->Write("JetPt_Area_Wide_10_40");
            JetPt_Area_Wide[7][0]->Write("JetPt_Area_Wide_40_80");

            Z_Area_Wide[0][0]->Write("Z_Area_Wide");
            Z_Area_Wide[1][0]->Write("Z_Area_Wide_0_10");
            Z_Area_Wide[2][0]->Write("Z_Area_Wide_10_20");
            Z_Area_Wide[3][0]->Write("Z_Area_Wide_20_30");
            Z_Area_Wide[4][0]->Write("Z_Area_Wide_30_40");
            Z_Area_Wide[6][0]->Write("Z_Area_Wide_10_40");
            Z_Area_Wide[7][0]->Write("Z_Area_Wide_40_80");

            ZPt_Area_Wide[0][0]->Write("ZPt_Area_Wide");
            ZPt_Area_Wide[1][0]->Write("ZPt_Area_Wide_0_10");
            ZPt_Area_Wide[2][0]->Write("ZPt_Area_Wide_10_20");
            ZPt_Area_Wide[3][0]->Write("ZPt_Area_Wide_20_30");
            ZPt_Area_Wide[4][0]->Write("ZPt_Area_Wide_30_40");
            ZPt_Area_Wide[6][0]->Write("ZPt_Area_Wide_10_40");
            ZPt_Area_Wide[7][0]->Write("ZPt_Area_Wide_40_80");

            JetPt_R_Area_Wide[0][0]->Write("JetPt_R_Area_Wide");
            JetPt_R_Area_Wide[1][0]->Write("JetPt_R_Area_Wide_0_10");
            JetPt_R_Area_Wide[2][0]->Write("JetPt_R_Area_Wide_10_20");
            JetPt_R_Area_Wide[3][0]->Write("JetPt_R_Area_Wide_20_30");
            JetPt_R_Area_Wide[4][0]->Write("JetPt_R_Area_Wide_30_40");
            JetPt_R_Area_Wide[6][0]->Write("JetPt_R_Area_Wide_10_40");
            JetPt_R_Area_Wide[7][0]->Write("JetPt_R_Area_Wide_40_80");

            JetPtAreaVsJetPtCS[0][0]->Write("JetPtAreaVsJetPtCS");
            JetPtAreaVsJetPtCS[1][0]->Write("JetPtAreaVsJetPtCS_0_10");
            JetPtAreaVsJetPtCS[2][0]->Write("JetPtAreaVsJetPtCS_10_20");
            JetPtAreaVsJetPtCS[3][0]->Write("JetPtAreaVsJetPtCS_20_30");
            JetPtAreaVsJetPtCS[4][0]->Write("JetPtAreaVsJetPtCS_30_40");
            JetPtAreaVsJetPtCS[6][0]->Write("JetPtAreaVsJetPtCS_10_40");
            JetPtAreaVsJetPtCS[7][0]->Write("JetPtAreaVsJetPtCS_40_80");

            JetAreaVsJetNConstArea[0][0]->Write("JetAreaVsJetNConstArea");
            JetAreaVsJetNConstArea[1][0]->Write("JetAreaVsJetNConstArea_0_10");
            JetAreaVsJetNConstArea[2][0]->Write("JetAreaVsJetNConstArea_10_20");
            JetAreaVsJetNConstArea[3][0]->Write("JetAreaVsJetNConstArea_20_30");
            JetAreaVsJetNConstArea[4][0]->Write("JetAreaVsJetNConstArea_30_40");
            JetAreaVsJetNConstArea[6][0]->Write("JetAreaVsJetNConstArea_10_40");
            JetAreaVsJetNConstArea[7][0]->Write("JetAreaVsJetNConstArea_40_80");

            JetAreaVsJetNConstCS[0][0]->Write("JetAreaVsJetNConstCS");
            JetAreaVsJetNConstCS[1][0]->Write("JetAreaVsJetNConstCS_0_10");
            JetAreaVsJetNConstCS[2][0]->Write("JetAreaVsJetNConstCS_10_20");
            JetAreaVsJetNConstCS[3][0]->Write("JetAreaVsJetNConstCS_20_30");
            JetAreaVsJetNConstCS[4][0]->Write("JetAreaVsJetNConstCS_30_40");
            JetAreaVsJetNConstCS[6][0]->Write("JetAreaVsJetNConstCS_10_40");
            JetAreaVsJetNConstCS[7][0]->Write("JetAreaVsJetNConstCS_40_80");

            JetPtAreaVsJetNConstArea[0][0]->Write("JetPtAreaVsJetNConstArea");
            JetPtAreaVsJetNConstArea[1][0]->Write("JetPtAreaVsJetNConstArea_0_10");
            JetPtAreaVsJetNConstArea[2][0]->Write("JetPtAreaVsJetNConstArea_10_20");
            JetPtAreaVsJetNConstArea[3][0]->Write("JetPtAreaVsJetNConstArea_20_30");
            JetPtAreaVsJetNConstArea[4][0]->Write("JetPtAreaVsJetNConstArea_30_40");
            JetPtAreaVsJetNConstArea[6][0]->Write("JetPtAreaVsJetNConstArea_10_40");
            JetPtAreaVsJetNConstArea[7][0]->Write("JetPtAreaVsJetNConstArea_40_80");

            JetPtCSVsJetNConstCS[0][0]->Write("JetPtCSVsJetNConstCS");
            JetPtCSVsJetNConstCS[1][0]->Write("JetPtCSVsJetNConstCS_0_10");
            JetPtCSVsJetNConstCS[2][0]->Write("JetPtCSVsJetNConstCS_10_20");
            JetPtCSVsJetNConstCS[3][0]->Write("JetPtCSVsJetNConstCS_20_30");
            JetPtCSVsJetNConstCS[4][0]->Write("JetPtCSVsJetNConstCS_30_40");
            JetPtCSVsJetNConstCS[6][0]->Write("JetPtCSVsJetNConstCS_10_40");
            JetPtCSVsJetNConstCS[7][0]->Write("JetPtCSVsJetNConstCS_40_80");

            JetCentVsJetNConst[0][0]->Write("JetCentVsJetNConst");

            JetRefMult[0][0]->Write("JetRefMult");
        }
    }
}

// void doSys()
// {
//     SYS = true;
//     fillSys(D0Sys_010, 1);
//     fillSys(D0Sys_1040, 2);
//     fillSys(D0Sys_4080, 3);
//     for (int i = 0; i < 166; i++)
//     {
//         cout << ">>> On systematic loop " << i + 1 << endl;
//         // sprintf(NAME,"systematics_new_ratios/Histograms3_D05GeV_%i.root",i);
//         cout << " >> Will save histograms as " << NAME << endl;
//         shuffleEff();
//         applyweights(1);
//     }
// }
