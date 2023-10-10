#include "../Config.h"

TFile *eff = new TFile("D0Eff.root","READ");
TH1F* hEff0 = (TH1F*) eff->Get("hEff0");
TH1F* hEff1 = (TH1F*) eff->Get("hEff1");
TH1F* hEff2 = (TH1F*) eff->Get("hEff2");
TH1F* hEff3 = (TH1F*) eff->Get("hEff3");
TH1F* hEff4 = (TH1F*) eff->Get("hEff4");

TH1F* hEff0_Sys = (TH1F*) hEff0->Clone("hEff0_Sys");
TH1F* hEff1_Sys = (TH1F*) hEff1->Clone("hEff1_Sys");
TH1F* hEff2_Sys = (TH1F*) hEff2->Clone("hEff2_Sys");
TH1F* hEff3_Sys = (TH1F*) hEff3->Clone("hEff3_Sys");
TH1F* hEff4_Sys = (TH1F*) hEff4->Clone("hEff4_Sys");

TFile *Eff_Sys = new TFile("Sys4Matt.root","READ");
TGraph *gSys_cent0_err0 = (TGraph*)Eff_Sys->Get("gSys_cent0_err0");
TGraph *gSys_cent0_err1 = (TGraph*)Eff_Sys->Get("gSys_cent0_err1");
TGraph *gSys_cent0_err2 = (TGraph*)Eff_Sys->Get("gSys_cent0_err2");
TGraph *gSys_cent0_err3 = (TGraph*)Eff_Sys->Get("gSys_cent0_err3");
TGraph *gSys_cent0_err4 = (TGraph*)Eff_Sys->Get("gSys_cent0_err4");
TGraph *gSys_cent0_err5 = (TGraph*)Eff_Sys->Get("gSys_cent0_err5");
TGraph *gSys_cent0_err6 = (TGraph*)Eff_Sys->Get("gSys_cent0_err6");
TGraph *gSys_cent0_err7 = (TGraph*)Eff_Sys->Get("gSys_cent0_err7");
TGraph *gSys_cent0_err8 = (TGraph*)Eff_Sys->Get("gSys_cent0_err8");
TGraph *gSys_cent0_err9 = (TGraph*)Eff_Sys->Get("gSys_cent0_err9");

TGraph *gSys_cent5_err0 = (TGraph*)Eff_Sys->Get("gSys_cent5_err0");
TGraph *gSys_cent5_err1 = (TGraph*)Eff_Sys->Get("gSys_cent5_err1");
TGraph *gSys_cent5_err2 = (TGraph*)Eff_Sys->Get("gSys_cent5_err2");
TGraph *gSys_cent5_err3 = (TGraph*)Eff_Sys->Get("gSys_cent5_err3");
TGraph *gSys_cent5_err4 = (TGraph*)Eff_Sys->Get("gSys_cent5_err4");
TGraph *gSys_cent5_err5 = (TGraph*)Eff_Sys->Get("gSys_cent5_err5");
TGraph *gSys_cent5_err6 = (TGraph*)Eff_Sys->Get("gSys_cent5_err6");
TGraph *gSys_cent5_err7 = (TGraph*)Eff_Sys->Get("gSys_cent5_err7");
TGraph *gSys_cent5_err8 = (TGraph*)Eff_Sys->Get("gSys_cent5_err8");
TGraph *gSys_cent5_err9 = (TGraph*)Eff_Sys->Get("gSys_cent5_err9");

TGraph *gSys_cent6_err0 = (TGraph*)Eff_Sys->Get("gSys_cent6_err0");
TGraph *gSys_cent6_err1 = (TGraph*)Eff_Sys->Get("gSys_cent6_err1");
TGraph *gSys_cent6_err2 = (TGraph*)Eff_Sys->Get("gSys_cent6_err2");
TGraph *gSys_cent6_err3 = (TGraph*)Eff_Sys->Get("gSys_cent6_err3");
TGraph *gSys_cent6_err4 = (TGraph*)Eff_Sys->Get("gSys_cent6_err4");
TGraph *gSys_cent6_err5 = (TGraph*)Eff_Sys->Get("gSys_cent6_err5");
TGraph *gSys_cent6_err6 = (TGraph*)Eff_Sys->Get("gSys_cent6_err6");
TGraph *gSys_cent6_err7 = (TGraph*)Eff_Sys->Get("gSys_cent6_err7");
TGraph *gSys_cent6_err8 = (TGraph*)Eff_Sys->Get("gSys_cent6_err8");
TGraph *gSys_cent6_err9 = (TGraph*)Eff_Sys->Get("gSys_cent6_err9");

TH1F *D0Sys_010 = new TH1F("D0Sys_010","D0Sys_010",PTBINS,edges);
TH1F *D0Sys_1040 = new TH1F("D0Sys_1040","D0Sys_1040",PTBINS,edges);
TH1F *D0Sys_4080 = new TH1F("D0Sys_480","D0Sys_4080",PTBINS,edges);

TFile *doub = new TFile("MisPID_SB_Final.root","READ");
TGraph *gDoub0 = (TGraph*)doub->Get("DoubleCounting_Cen_0_SB");
TGraph *gDoub1 = (TGraph*)doub->Get("DoubleCounting_Cen_1_SB");
TGraph *gDoub2 = (TGraph*)doub->Get("DoubleCounting_Cen_2_SB");
TGraph *gDoub3 = (TGraph*)doub->Get("DoubleCounting_Cen_3_SB");
bool useRefmult = true;
bool useEff = true; // aPply eff corr if true
bool useDouble = true; //Apply double corr. if true
// bool usePt = true; // Apply Jet pT 5 3 geV if true
bool usePt = true; // Apply Jet pT > 3 geV if true
bool useHighD0Pt = true; // Appl D0 pt > 3 GeV if true
bool useLowD0Pt = false;// for D0 1-3 GeV //changing it to D0 1-10 GeV now 
bool useD0Eta = false; //D0 eta in 0.6
bool useJetArea = false; //Jet area cut
bool useBKG = false;//FOr BKG distributions
bool useConst = false; // (D0-Jet) Pt <0 thronw away if true
bool useZ = false;// Cut events with D0 z <0 or > 1
bool SYS = false; //Run systematic uncertainties
char NAME[100] = "Histograms3_D05_10GeV_NewBinning.root";
int SAVE = 1;

double getEff(int mCen,float mPt){
    double eff=0;

    // Note I am using the efficiency histograms with "_Sys" here. If running applyWeights() only, this doesn't matter
    // If running doSys(), this will shuffle the efficinecy with a Gaussian each iteration in the _Sys histograms.
    // To not have to change the code here, _Sys is used all the time but is equivelent to original histograms in nominal running
    if(mCen>=0 && mCen<10)eff=hEff0_Sys->GetBinContent(hEff0_Sys->FindBin(mPt));
    else if(mCen>=10 && mCen<20)eff=hEff1_Sys->GetBinContent(hEff1_Sys->FindBin(mPt));
    else if(mCen>=20 && mCen<40)eff=hEff2_Sys->GetBinContent(hEff2_Sys->FindBin(mPt));
    else if(mCen>=40 && mCen<60)eff=hEff3_Sys->GetBinContent(hEff3_Sys->FindBin(mPt));
    else if(mCen>=60 && mCen<80)eff=(2./3.)*hEff4_Sys->GetBinContent(hEff4_Sys->FindBin(mPt));// efficiency is from paper which has scale factor
    return eff;
}
double getDoubleCount(int mCen,float mPt){
    double dcount=1;
    if(mCen ==0 && mCen == 80)dcount = gDoub0->Eval(mPt);
    else if(mCen>=0 && mCen<10)dcount = gDoub1->Eval(mPt);
    else if(mCen>=10 && mCen<40)dcount = gDoub2->Eval(mPt);
    else if(mCen>=40 && mCen<80)dcount = gDoub3->Eval(mPt);
    return dcount;
}

void DN(TH1F *R, int col, int mar, double sum){
    for(int i = 1;i<R->GetNbinsX()+1;i++){
	double val = R->GetBinContent(i);
	double er = R->GetBinError(i);
	double width = R->GetBinWidth(i);
	R->SetBinContent(i,val/width/sum);
	R->SetBinError(i,er/width/sum);
    }
//    R->Scale(1./sum);
    R->SetLineColor(col);
    R->SetMarkerColor(col);
    R->SetMarkerStyle(mar);
}

double fillHistos(char file[100], TH1F* D0Pt,TH1F* JetPt, TH1F* Z, TH1F* R, TH2F* ZPt,int cen_low,int cen_high, TH1F* hArea, TH1F* JetEta, TH1F* D0Eta, TH1F* JetPhi, TH1F* D0Phi, TH2F* D0JetEta, TH1F *Eta, TH2F*JetPt_R){
    if(!SYS)cout <<" >>> On fill with: " << file << endl;
    int sumw=0;
    TFile *f_D = new TFile(file);
    TChain* ch = (TChain*) f_D->Get("Signal_sw");
    float mM; float mJA;float weight; float mPt; float mJPt; float mR; float mZ; int mCen;float mW;float mEta;float mJEta;float mPhi; float mJPhi;
    ch->SetBranchAddress( "mM" , &mM );
    ch->SetBranchAddress( "sWeight" , &weight );
    ch->SetBranchAddress( "mPt" , &mPt );
    ch->SetBranchAddress( "mJPt" , &mJPt );
    ch->SetBranchAddress( "mR" , &mR );
    ch->SetBranchAddress( "mZ" , &mZ );
    ch->SetBranchAddress( "mCen" , &mCen );
    ch->SetBranchAddress( "mWeight" , &mW );
    ch->SetBranchAddress( "mJetArea" , &mJA );
    ch->SetBranchAddress( "mEta" , &mEta );
    ch->SetBranchAddress( "mJEta" , &mJEta );

    double sigma=0;
    if(cen_low==0 && cen_high==10)sigma= 20;//0.156415;
    if(cen_low==10 && cen_high==20)sigma= 16;//0.143552;
    if(cen_low==20 && cen_high==30)sigma= 16;//0.135034;
    if(cen_low==30 && cen_high==40)sigma= 16;//0.120011;
    if(cen_low==40 && cen_high==60)sigma= 10;//0.100353;
    if(cen_low==10 && cen_high==40)sigma= 16;//00.136703;
    if(cen_low==40 && cen_high==80)sigma= 10;//00.0972255;
    int loop = ch->GetEntries();
    
    
    for(int i =0;i<loop;i++){
        ch->GetEntry(i);
        // if (mJPt < 0) cout << mJPt << "\t" << mPt << endl;
        double eff = getEff(mCen,mPt);
	double dcount = getDoubleCount(mCen,mPt);
        if(mCen<cen_low || mCen>=cen_high)continue;
        if(useHighD0Pt && mPt<=5)continue;
	if(useLowD0Pt && (mPt < 1 || mPt>10))continue;
	if(usePt && mJPt<0)continue;

    // if (mZ > 1.) cout << mJPt << "\t" << mPt << "\t" << mZ << endl;

        if(useD0Eta && fabs(mJEta)>0.3)continue;
	if(useJetArea && mJA < 0.35)continue;
	if(useConst && (mPt - mJPt)>0)continue;
	if(useZ && (mZ <0. || mZ > 1.))continue;

    // if (mZ > 1.) cout << " ======== " << mJPt << "\t" << mPt << "\t" << mZ << endl;
	if(useEff && useDouble) weight *= (1.-dcount) / eff ;
        else if(useDouble && !useEff)  weight *= (1.-dcount);
        else if(!useDouble && useEff)  weight *= 1 / eff;
        if(useRefmult)weight *= mW;
        if(useBKG) weight=1.-weight;
        sumw+=weight;
        D0Pt->Fill(mPt,weight);
        JetPt->Fill(mJPt,weight);
        Z->Fill(mZ,weight);
	hArea->Fill(mJA,weight);
        R->Fill(mR,weight);
        ZPt->Fill(mJPt,mZ,weight);
	JetEta->Fill(mJEta,weight);
	D0Eta->Fill(mEta,weight);
        JetPhi->Fill(mJPhi,weight);
        D0Phi->Fill(mPhi,weight);
	D0JetEta->Fill(mEta,mJEta,weight);
	Eta->Fill(mJEta-mEta);
	JetPt_R->Fill(mJPt,mR,weight);
	
    }
    f_D->Close();
    return sumw;
}
double fillHistos2(char file[100], TH1F* D0Pt,TH1F* JetPt, TH1F* Z, TH1F* R, TH2F* ZPt,TH2F* D0JPt,TH2F *hard,int cen_low,int cen_high, TH1F* hArea, TH1F* JetEta, TH1F* D0Eta, TH1F* JetPhi, TH1F* D0Phi, TH2F* D0JetEta, TH1F *Eta, TH2F*JetPt_R){
    if(!SYS)cout <<" >>> On fill with: " << file << endl;
    int sumw=0;
    TFile *f_D = new TFile(file);
    TChain* ch = (TChain*) f_D->Get("Signal_sw");
    float mM; float mJA;float weight; float mPt; float mJPt; float mR; float mZ; int mCen;float mHPt;float mW;float mEta;float mJEta;float mPhi; float mJPhi;
    ch->SetBranchAddress( "mM" , &mM );
    ch->SetBranchAddress( "sWeight" , &weight );
    ch->SetBranchAddress( "mPt" , &mPt );
    ch->SetBranchAddress( "mJPt" , &mJPt );
    ch->SetBranchAddress( "mR" , &mR );
    ch->SetBranchAddress( "mZ" , &mZ );
    ch->SetBranchAddress( "mCen" , &mCen );
    ch->SetBranchAddress( "mHPt" , &mHPt );
    ch->SetBranchAddress( "mWeight" , &mW );
    ch->SetBranchAddress( "mJetArea" , &mJA );
    ch->SetBranchAddress( "mEta" , &mEta );
    ch->SetBranchAddress( "mJEta" , &mJEta );
    double sigma=0;
    if(cen_low==0 && cen_high==10)sigma= 20;//0.156415;                                                                                                                                                                                            
    if(cen_low==10 && cen_high==20)sigma= 16;//0.143552;                                                                                                                                                                                           
    if(cen_low==20 && cen_high==30)sigma= 16;//0.135034;                                                                                                                                                                                           
    if(cen_low==30 && cen_high==40)sigma= 16;//0.120011;                                                                                                                                                                                           
    if(cen_low==40 && cen_high==60)sigma= 10;//0.100353;                                                                                                                                                                                           
    if(cen_low==10 && cen_high==40)sigma= 16;//00.136703;                                                                                                                                                                                          
    if(cen_low==40 && cen_high==80)sigma= 10;//00.0972255;           
    int loop = ch->GetEntries();
    for(int i =0;i<loop;i++){
        ch->GetEntry(i);
        double eff = getEff(mCen,mPt);
	double dcount = getDoubleCount(mCen,mPt);
        if(mCen<cen_low || mCen>=cen_high)continue;
	if(useHighD0Pt && mPt<5)continue;
	if(useLowD0Pt && (mPt < 1 || mPt>10))continue;
        if(usePt && mJPt<0)continue;
	if(useD0Eta && fabs(mJEta)>0.6)continue;
	if(useJetArea && mJA < 0.35)continue;
	if(useConst && (mPt - mJPt)>0)continue;
	if(useZ&& (mZ <0. || mZ > 1.))continue;
	if(useEff && useDouble) weight *= (1.-dcount) / eff ;
        else if(useDouble && !useEff)  weight *= (1.-dcount);
        else if(!useDouble && useEff)  weight *= 1 / eff;
        if(useRefmult)weight *= mW;
        if(useBKG) weight=1.-weight;
        sumw+=weight;
        D0Pt->Fill(mPt,weight);
        D0JPt->Fill(mJPt,mPt,weight);
        JetPt->Fill(mJPt,weight);
        Z->Fill(mZ,weight);
        R->Fill(mR,weight);
        ZPt->Fill(mJPt,mZ,weight);
        hard->Fill(mJPt,mHPt/mPt,weight);
	JetEta->Fill(mJEta,weight);
	hArea->Fill(mJA,weight);
	D0Eta->Fill(mEta,weight);
        JetPhi->Fill(mJPhi,weight);
        D0Phi->Fill(mPhi,weight);
	D0JetEta->Fill(mEta,mJEta,weight);
	Eta->Fill(mJEta-mEta);
	JetPt_R->Fill(mJPt,mR,weight);

    }
    f_D->Close();
    return sumw;
}

double fillHistos1(char file[100],TH2F* PT2D, TH1F* D0Pt,TH1F* JetPt, TH1F* Z, TH1F* R,TH2F* ZPt, int cen_low,int cen_high, TH1F* centrality, TH1F* hArea, TH1F* JetEta, TH1F* D0Eta, TH1F* JetPhi, TH1F* D0Phi, TH2F* D0JetEta, TH1F *Eta, TH2F*JetPt_R){
    if(!SYS)cout <<" >>> On fill with: " << file << endl;
    int sumw=0;
    TFile *f_D = new TFile(file);
    TChain* ch = (TChain*) f_D->Get("Signal_sw");
    float mM; float mJA; float weight; float mPt; float mJPt; float mR; float mZ; int mCen;float mW;float mEta; float mJEta;float mPhi; float mJPhi;
    ch->SetBranchAddress( "mM" , &mM );
    ch->SetBranchAddress( "sWeight" , &weight );
    ch->SetBranchAddress( "mPt" , &mPt );
    ch->SetBranchAddress( "mJPt" , &mJPt );
    ch->SetBranchAddress( "mR" , &mR );
    ch->SetBranchAddress( "mZ" , &mZ );
    ch->SetBranchAddress( "mCen" , &mCen );
    ch->SetBranchAddress( "mWeight" , &mW );
    ch->SetBranchAddress( "mJetArea" , &mJA );
    ch->SetBranchAddress( "mEta" , &mEta );
    ch->SetBranchAddress( "mJEta" , &mJEta );    
    double sigma=0;
    if(cen_low==0 && cen_high==10)sigma= 20;//0.156415;                                                                                                                                                                                            
    if(cen_low==10 && cen_high==20)sigma= 16;//0.143552;                                                                                                                                                                                           
    if(cen_low==20 && cen_high==30)sigma= 16;//0.135034;                                                                                                                                                                                           
    if(cen_low==30 && cen_high==40)sigma= 16;//0.120011;                                                                                                                                                                                           
    if(cen_low==40 && cen_high==60)sigma= 10;//0.100353;                                                                                                                                                                                           
    if(cen_low==10 && cen_high==40)sigma= 16;//00.136703;                                                                                                                                                                                          
    if(cen_low==40 && cen_high==80)sigma= 10;//00.0972255;           

    int loop = ch->GetEntries();
    for(int i =0;i<loop;i++){
        ch->GetEntry(i);
        double eff = getEff(mCen,mPt);
 	double dcount = getDoubleCount(mCen,mPt);
	if(useHighD0Pt && mPt<5)continue;
	if(useLowD0Pt && (mPt < 1 || mPt>10))continue;
        if(usePt && mJPt<0)continue;

    if (mZ > 1.) cout << mJPt << "\t" << mPt << "\t" << mZ << endl;


	if(useD0Eta && fabs(mJEta)>0.6)continue;
	if(useJetArea && mJA < 0.35)continue;
	if(useConst && (mPt - mJPt)>0)continue;
	// if(useZ&& (mZ <0. || mZ > 1.))continue;

    if (mZ > 1.) cout << " ======== " << mJPt << "\t" << mPt << "\t" << mZ << endl;
	if(useEff && useDouble) weight *= (1.-dcount) / eff ;
        else if(useDouble && !useEff)  weight *= (1.-dcount);
        else if(!useDouble && useEff)  weight *= 1 / eff;
        if(useRefmult)weight *= mW;
        if(useBKG) weight=1.-weight;
        sumw+=weight;


	D0Pt->Fill(mPt,weight);
        JetPt->Fill(mJPt,weight);
        Z->Fill(mZ,weight);
        R->Fill(mR,weight);
        ZPt->Fill(mJPt,mZ,weight);
        PT2D->Fill(mJPt,mPt,weight);
        centrality->Fill(mCen+0.1,weight);
	hArea->Fill(mJA,weight);
	JetEta->Fill(mJEta,weight);
	D0Eta->Fill(mEta,weight);
	JetPhi->Fill(mJPhi,weight);
        D0Phi->Fill(mPhi,weight);
	D0JetEta->Fill(mEta,mJEta,weight);
	Eta->Fill(mJEta-mEta);
	JetPt_R->Fill(mJPt,mR,weight);
}
    f_D->Close();
    return sumw;
}
void shuffleEff(){
    for(int i = 1; i <hEff0->GetNbinsX()+1;i++){
        double val0 = hEff0->GetBinContent(i);
        double val1 = hEff1->GetBinContent(i);
        double val2 = hEff2->GetBinContent(i);
        double val3 = hEff3->GetBinContent(i);
        double val4 = hEff4->GetBinContent(i);
        double err0 = D0Sys_010->GetBinContent(i);//hEff0->GetBinError(i);
        double err1 = D0Sys_1040->GetBinContent(i);//hEff1->GetBinError(i);
        double err2 = D0Sys_1040->GetBinContent(i);//hEff2->GetBinError(i);
        double err3 = D0Sys_1040->GetBinContent(i);//hEff3->GetBinError(i);
        double err4 = D0Sys_4080->GetBinContent(i);//hEff4->GetBinError(i);
        float shift0 = gRandom->Gaus(0,err0);
        float shift1 = gRandom->Gaus(0,err1);
        float shift2 = gRandom->Gaus(0,err2);
        float shift3 = gRandom->Gaus(0,err3);
        float shift4 = gRandom->Gaus(0,err4);
	hEff0_Sys->SetBinContent(i,val0/(1+shift0));
        hEff1_Sys->SetBinContent(i,val1/(1+shift1));
        hEff2_Sys->SetBinContent(i,val2/(1+shift2));
        hEff3_Sys->SetBinContent(i,val3/(1+shift3));
        hEff4_Sys->SetBinContent(i,val4/(1+shift4));
    }
}

void fillSys(TH1* h, int cen){
    
    for(int i = 1;i<h->GetNbinsX()+1;i++){
	double v1,v2,v3,v4,v5,v6,v7,v8,v9,v10;
	if(cen==1){
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
	    
	}else if(cen==2){
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
	}else if(cen==3){

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
	//Ignore v3, which comes from yield extraction
	//double total = sqrt(v1*v1+v2*v2+v4*v4+v5*v5+v6*v6+v7*v7+v8*v8+v9*v9+v10*v10);
	//for ratios take out correalted ones
	double total = sqrt(v1*v1+v2*v2+v4*v4+v5*v5+v6*v6+v7*v7+v10*v10);

	h->SetBinContent(i,total);
	cout<< "Cen " << cen << " bin " << i << " total "<< total << endl;
    
    }
}

void applyWeightsNew(int save = SAVE){
    if(!SYS)gROOT->ProcessLine(".x ../myStyle.C");
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    double sumw=0;
    double sumw1=0;
    double sumw2=0;
    double sumw3=0;  
    double sumw4=0;
    double sumw5=0;
    double sumw6=0;
    double sumw7=0;
    double sumw8=0;
    double sumw9=0;
    double sumw10=0;
    double binning[6] = {0,0.05,0.1,0.2,0.3,0.4};
    double binning1[9] = {5,7,9,11,13,15,20,30,40};
    double binning2[11] = {0,1,2,3,4,5,6,8,10};
    TH2F *PT2D = new TH2F("PT2D","PT2D",35,5,40,50,0,10);//8,binning1,10,binning2);
    TH1F *hArea = new TH1F("hArea","hArea",100,0,1);
    TH1F *D0Pt = new TH1F("D0Pt","D0Pt",PTBINS,edges);
    TH1F *JetPt = new TH1F("JetPt","JetPt",nbins_jpt,binning_jpt);
    TH2F *JetPt_R = new TH2F("JetPt_R","JetPt_R",nBinsJetPt_R,JetPtBins_R,nBinsdR,dRBins);
    TH1F *Z = new TH1F("Z","Z",nbins_z,binning_z);
    TH2F *ZPt = new TH2F("ZPt","ZPt",nbins_jpt,binning_jpt,nbins_z,binning_z);
    TH1F *R = new TH1F("R","R",nBinsdR,dRBins);
    TH1F *D0Pt_1 = new TH1F("D0Pt_1","D0Pt_1",PTBINS,edges);//50,0,50);
    TH1F *JetPt_1 = new TH1F("JetPt_1","JetPt_1",nbins_jpt,binning_jpt);
    TH2F *JetPt_R_1 = new TH2F("JetPt_R_1","JetPt_R_1",nBinsJetPt_R,JetPtBins_R,nBinsdR,dRBins);
    TH1F *Z_1 = new TH1F("Z_1","Z_1",nbins_z,binning_z);
    TH1F *R_1 = new TH1F("R_1","R_1",nBinsdR,dRBins);
    TH2F *ZPt_1 = new TH2F("ZPt_1","ZPt_1",nbins_jpt,binning_jpt,nbins_z,binning_z);
    TH1F *D0Pt_2 = new TH1F("D0Pt_2","D0Pt_2",PTBINS,edges);//50,0,50);
    TH1F *JetPt_2 = new TH1F("JetPt_2","JetPt_2",nbins_jpt,binning_jpt);
    TH2F *JetPt_R_2 = new TH2F("JetPt_R_2","JetPt_R_2",nBinsJetPt_R,JetPtBins_R,nBinsdR,dRBins);
    TH1F *Z_2 = new TH1F("Z_2","Z_2",nbins_z,binning_z);
    TH1F *R_2 = new TH1F("R_2","R_2",nBinsdR,dRBins);
    TH2F *ZPt_2 = new TH2F("ZPt_2","ZPt_2",nbins_jpt,binning_jpt,nbins_z,binning_z);
    TH1F *D0Pt_3 = new TH1F("D0Pt_3","D0Pt_3",PTBINS,edges);//50,0,50);
    TH1F *JetPt_3 = new TH1F("JetPt_3","JetPt_3",nbins_jpt,binning_jpt);
    TH2F *JetPt_R_3 = new TH2F("JetPt_R_3","JetPt_R_3",nBinsJetPt_R,JetPtBins_R,nBinsdR,dRBins);
    TH1F *Z_3 = new TH1F("Z_3","Z_3",nbins_z,binning_z);
  TH1F *R_3 = new TH1F("R_3","R_3",nBinsdR,dRBins);
  TH2F *ZPt_3 = new TH2F("ZPt_3","ZPt_3",nbins_jpt,binning_jpt,nbins_z,binning_z);
  TH1F *D0Pt_4 = new TH1F("D0Pt_4","D0Pt_4",PTBINS,edges);//50,0,50);
  TH1F *JetPt_4 = new TH1F("JetPt_4","JetPt_4",nbins_jpt,binning_jpt);
  TH2F *JetPt_R_4 = new TH2F("JetPt_R_4","JetPt_R_4",nBinsJetPt_R,JetPtBins_R,nBinsdR,dRBins);
  TH1F *Z_4 = new TH1F("Z_4","Z_4",nbins_z,binning_z);
  TH1F *R_4 = new TH1F("R_4","R_4",nBinsdR,dRBins);
  TH2F *ZPt_4 = new TH2F("ZPt_4","ZPt_4",nbins_jpt,binning_jpt,nbins_z,binning_z);
  TH1F *D0Pt_5 = new TH1F("D0Pt_5","D0Pt_5",PTBINS,edges);//50,0,50);
  TH2F *D0JPt_5 = new TH2F("D0JPt_5","D0JPt_5",nbins_jpt,binning_jpt,PTBINS,edges);//50,0,50);   
  TH1F *JetPt_5 = new TH1F("JetPt_5","JetPt_5",nbins_jpt,binning_jpt);
  TH2F *JetPt_R_5 = new TH2F("JetPt_R_5","JetPt_R_5",nBinsJetPt_R,JetPtBins_R,nBinsdR,dRBins);
  TH1F *Z_5 = new TH1F("Z_5","Z_5",nbins_z,binning_z);
  TH1F *R_5 = new TH1F("R_5","R_5",nBinsdR,dRBins);
  TH2F *ZPt_5 = new TH2F("ZPt_5","ZPt_5",nbins_jpt,binning_jpt,nbins_z,binning_z);
  TH2F *Hard2D2_5 = new TH2F("r_Hard2D2_5","r_Hard2D2_5",100,0,50,100,0,5);
  TH1F *D0Pt_6 = new TH1F("D0Pt_6","D0Pt_6",PTBINS,edges);//50,0,50);                                                                                                                                                          
  TH1F *JetPt_6 = new TH1F("JetPt_6","JetPt_6",nbins_jpt,binning_jpt);
  TH2F *JetPt_R_6 = new TH2F("JetPt_R_6","JetPt_R_6",nBinsJetPt_R,JetPtBins_R,nBinsdR,dRBins);
  TH1F *Z_6 = new TH1F("Z_6","Z_6",nbins_z,binning_z);
  TH1F *R_6 = new TH1F("R_6","R_6",nBinsdR,dRBins);
  TH2F *ZPt_6 = new TH2F("ZPt_6","ZPt_6",nbins_jpt,binning_jpt,nbins_z,binning_z);

  TH1F *D0Pt_7 = new TH1F("D0Pt_7","D0Pt_7",PTBINS,edges);//50,0,50);                                                                                                                                                          
  TH1F *JetPt_7 = new TH1F("JetPt_7","JetPt_7",nbins_jpt,binning_jpt);
  TH1F *Z_7 = new TH1F("Z_7","Z_7",nbins_z,binning_z);
  TH1F *R_7 = new TH1F("R_7","R_7",nBinsdR,dRBins);
  TH2F *ZPt_7 = new TH2F("ZPt_7","ZPt_7",nbins_jpt,binning_jpt,nbins_z,binning_z);  
  TH2F *JetPt_R_7 = new TH2F("JetPt_R_7","JetPt_R_7",nBinsJetPt_R,JetPtBins_R,nBinsdR,dRBins);

  TH1F *D0Pt_8 = new TH1F("D0Pt_8","D0Pt_8",PTBINS,edges);//50,0,50);             
  TH1F *JetPt_8 = new TH1F("JetPt_8","JetPt_8",njpt_bins,jetpt_low, jetpt_high);
  TH1F *Z_8 = new TH1F("Z_8","Z_8",nz_bins,z_low, z_high);
  TH1F *R_8 = new TH1F("R_8","R_8",nBinsdR,dRBins);
  TH2F *ZPt_8 = new TH2F("ZPt_8","ZPt_8",njpt_bins,jetpt_low, jetpt_high,nz_bins,z_low, z_high);  
  TH2F *JetPt_R_8 = new TH2F("JetPt_R_8","JetPt_R_8",njpt_bins,jetpt_low, jetpt_high,nBinsdR,dRBins);

  TH1F *D0Pt_9 = new TH1F("D0Pt_9","D0Pt_9",PTBINS,edges);//50,0,50);             
  TH1F *JetPt_9 = new TH1F("JetPt_9","JetPt_9",njpt_bins,jetpt_low, jetpt_high);
  TH1F *Z_9 = new TH1F("Z_9","Z_9",nz_bins,z_low, z_high);
  TH1F *R_9 = new TH1F("R_9","R_9",nBinsdR,dRBins);
  TH2F *ZPt_9 = new TH2F("ZPt_9","ZPt_9",njpt_bins,jetpt_low, jetpt_high,nz_bins,z_low, z_high);  
  TH2F *JetPt_R_9 = new TH2F("JetPt_R_9","JetPt_R_9",njpt_bins,jetpt_low, jetpt_high,nBinsdR,dRBins);

  TH1F *D0Pt_10 = new TH1F("D0Pt_10","D0Pt_10",PTBINS,edges);//50,0,50);             
  TH1F *JetPt_10 = new TH1F("JetPt_10","JetPt_10",njpt_bins,jetpt_low, jetpt_high);
  TH1F *Z_10 = new TH1F("Z_10","Z_10",nz_bins,z_low, z_high);
  TH1F *R_10 = new TH1F("R_10","R_10",nBinsdR,dRBins);
  TH2F *ZPt_10 = new TH2F("ZPt_10","ZPt_10",njpt_bins,jetpt_low, jetpt_high,nz_bins,z_low, z_high);  
  TH2F *JetPt_R_10 = new TH2F("JetPt_R_10","JetPt_R_10",njpt_bins,jetpt_low, jetpt_high,nBinsdR,dRBins);



  TH1F *Eta = new TH1F("Eta","D0Eta",200,-1,1);
  TH1F *Eta_1 = new TH1F("Eta_1","Eta_1",200,-1,1);
  TH1F *Eta_2 = new TH1F("Eta_2","Eta_2",200,-1,1);
  TH1F *Eta_3 = new TH1F("Eta_3","Eta_3",200,-1,1);
  TH1F *Eta_4 = new TH1F("Eta_4","Eta_4",200,-1,1);
  TH1F *Eta_5 = new TH1F("Eta_5","Eta_5",200,-1,1);
  TH1F *Eta_6 = new TH1F("Eta_6","Eta_6",200,-1,1);
  TH1F *Eta_7 = new TH1F("Eta_7","Eta_7",200,-1,1);
  TH1F *Eta_8 = new TH1F("Eta_8","Eta_8",200,-1,1);
  TH1F *Eta_9 = new TH1F("Eta_9","Eta_8",200,-1,1);
  TH1F *Eta_10 = new TH1F("Eta_10","Eta_10",200,-1,1);

  TH1F *D0Eta = new TH1F("D0Eta","D0Eta",200,-1,1); 
  TH1F *D0Eta_1 = new TH1F("D0Eta_1","D0Eta_1",200,-1,1);
  TH1F *D0Eta_2 = new TH1F("D0Eta_2","D0Eta_2",200,-1,1);
  TH1F *D0Eta_3 = new TH1F("D0Eta_3","D0Eta_3",200,-1,1);
  TH1F *D0Eta_4 = new TH1F("D0Eta_4","D0Eta_4",200,-1,1);
  TH1F *D0Eta_5 = new TH1F("D0Eta_5","D0Eta_5",200,-1,1);
  TH1F *D0Eta_6 = new TH1F("D0Eta_6","D0Eta_6",200,-1,1);
  TH1F *D0Eta_7 = new TH1F("D0Eta_7","D0Eta_7",200,-1,1);
  TH1F *D0Eta_8 = new TH1F("D0Eta_8","D0Eta_8",200,-1,1);
  TH1F *D0Eta_9 = new TH1F("D0Eta_9","D0Eta_9",200,-1,1);
  TH1F *D0Eta_10 = new TH1F("D0Eta_10","D0Eta_10",200,-1,1);


  TH1F *JetEta = new TH1F("JetEta","JetEta",200,-1,1);
  TH1F *JetEta_1 = new TH1F("JetEta_1","JetEta_1",200,-1,1);
  TH1F *JetEta_2 = new TH1F("JetEta_2","JetEta_2",200,-1,1);
  TH1F *JetEta_3 = new TH1F("JetEta_3","JetEta_3",200,-1,1);
  TH1F *JetEta_4 = new TH1F("JetEta_4","JetEta_4",200,-1,1);
  TH1F *JetEta_5 = new TH1F("JetEta_5","JetEta_5",200,-1,1);
  TH1F *JetEta_6 = new TH1F("JetEta_6","JetEta_6",200,-1,1);
  TH1F *JetEta_7 = new TH1F("JetEta_7","JetEta_7",200,-1,1);
  TH1F *JetEta_8 = new TH1F("JetEta_8","JetEta_8",200,-1,1);
  TH1F *JetEta_9 = new TH1F("JetEta_9","JetEta_9",200,-1,1);
  TH1F *JetEta_10 = new TH1F("JetEta_10","JetEta_10",200,-1,1);

  TH2F *D0JetEta = new TH2F("D0JetEta","D0JetEta",200,-1,1,200,-1,1);
  TH2F *D0JetEta_1 = new TH2F("D0JetEta_1","D0JetEta_1",200,-1,1,200,-1,1);
  TH2F *D0JetEta_2 = new TH2F("D0JetEta_2","D0JetEta_2",200,-1,1,200,-1,1);
  TH2F *D0JetEta_3 = new TH2F("D0JetEta_3","D0JetEta_3",200,-1,1,200,-1,1);
  TH2F *D0JetEta_4 = new TH2F("D0JetEta_4","D0JetEta_4",200,-1,1,200,-1,1);
  TH2F *D0JetEta_5 = new TH2F("D0JetEta_5","D0JetEta_5",200,-1,1,200,-1,1);
  TH2F *D0JetEta_6 = new TH2F("D0JetEta_6","D0JetEta_6",200,-1,1,200,-1,1);
  TH2F *D0JetEta_7 = new TH2F("D0JetEta_7","D0JetEta_7",200,-1,1,200,-1,1);
  TH2F *D0JetEta_8 = new TH2F("D0JetEta_8","D0JetEta_8",200,-1,1,200,-1,1);
  TH2F *D0JetEta_9 = new TH2F("D0JetEta_9","D0JetEta_9",200,-1,1,200,-1,1);
  TH2F *D0JetEta_10 = new TH2F("D0JetEta_10","D0JetEta_10",200,-1,1,200,-1,1);


  TH1F *D0Phi = new TH1F("D0Phi","D0Phi",200,-1,1);
  TH1F *D0Phi_1 = new TH1F("D0Phi_1","D0Phi_1",200,-1,1);
  TH1F *D0Phi_2 = new TH1F("D0Phi_2","D0Phi_2",200,-1,1);
  TH1F *D0Phi_3 = new TH1F("D0Phi_3","D0Phi_3",200,-1,1);
  TH1F *D0Phi_4 = new TH1F("D0Phi_4","D0Phi_4",200,-1,1);
  TH1F *D0Phi_5 = new TH1F("D0Phi_5","D0Phi_5",200,-1,1);
  TH1F *D0Phi_6 = new TH1F("D0Phi_6","D0Phi_6",200,-1,1);
  TH1F *D0Phi_7 = new TH1F("D0Phi_7","D0Phi_7",200,-1,1);
  TH1F *D0Phi_8 = new TH1F("D0Phi_8","D0Phi_8",200,-1,1);
  TH1F *D0Phi_9 = new TH1F("D0Phi_9","D0Phi_9",200,-1,1);
  TH1F *D0Phi_10 = new TH1F("D0Phi_10","D0Phi_10",200,-1,1);

  TH1F *JetPhi = new TH1F("JetPhi","JetPhi",200,-1,1);
  TH1F *JetPhi_1 = new TH1F("JetPhi_1","JetPhi_1",200,-1,1);
  TH1F *JetPhi_2 = new TH1F("JetPhi_2","JetPhi_2",200,-1,1);
  TH1F *JetPhi_3 = new TH1F("JetPhi_3","JetPhi_3",200,-1,1);
  TH1F *JetPhi_4 = new TH1F("JetPhi_4","JetPhi_4",200,-1,1);
  TH1F *JetPhi_5 = new TH1F("JetPhi_5","JetPhi_5",200,-1,1);
  TH1F *JetPhi_6 = new TH1F("JetPhi_6","JetPhi_6",200,-1,1);
  TH1F *JetPhi_7 = new TH1F("JetPhi_7","JetPhi_7",200,-1,1);
  TH1F *JetPhi_8 = new TH1F("JetPhi_8","JetPhi_8",200,-1,1);
  TH1F *JetPhi_9 = new TH1F("JetPhi_9","JetPhi_9",200,-1,1);
  TH1F *JetPhi_10 = new TH1F("JetPhi_10","JetPhi_10",200,-1,1);



  TH1F *hArea_1 = new TH1F("hArea_1","hArea_1",100,0,1);
  TH1F *hArea_2 = new TH1F("hArea_2","hArea_2",100,0,1);
  TH1F *hArea_3 = new TH1F("hArea_3","hArea_3",100,0,1);
  TH1F *hArea_4 = new TH1F("hArea_4","hArea_4",100,0,1);
  TH1F *hArea_5 = new TH1F("hArea_5","hArea_5",100,0,1);
  TH1F *hArea_6 = new TH1F("hArea_6","hArea_6",100,0,1);
  TH1F *hArea_7 = new TH1F("hArea_7","hArea_7",100,0,1);
  TH1F *hArea_8 = new TH1F("hArea_8","hArea_8",100,0,1);
  TH1F *hArea_9 = new TH1F("hArea_9","hArea_9",100,0,1);
  TH1F *hArea_10 = new TH1F("hArea_10","hArea_10",100,0,1);


  TH1F *centrality = new TH1F("centrality","centrality",20,0,100);

  // sumw+=fillHistos1("../sPlot/sWeight_PT_1_3_Cen_0_80.root",PT2D,D0Pt,JetPt,Z,R,ZPt,0,80,centrality,hArea,JetEta,D0Eta,JetPhi,D0Phi,D0JetEta,Eta,JetPt_R);
  // sumw+=fillHistos1("../sPlot/sWeight_PT_3_10_Cen_0_80.root",PT2D,D0Pt,JetPt,Z,R,ZPt,0,80,centrality,hArea,JetEta,D0Eta,JetPhi,D0Phi,D0JetEta,Eta,JetPt_R);
  sumw+=fillHistos1("../sPlot/sWeight_PT_1_10_Cen_0_80.root",PT2D,D0Pt,JetPt,Z,R,ZPt,0,80,centrality,hArea,JetEta,D0Eta,JetPhi,D0Phi,D0JetEta,Eta,JetPt_R);

  // sumw1+=fillHistos("../sPlot/sWeight_PT_1_3_Cen_0_10.root",D0Pt_1,JetPt_1,Z_1,R_1,ZPt_1,0,10,hArea_1,JetEta_1,D0Eta_1,JetPhi_1,D0Phi_1,D0JetEta_1,Eta_1,JetPt_R_1);
  // sumw1+=fillHistos("../sPlot/sWeight_PT_3_10_Cen_0_10.root",D0Pt_1,JetPt_1,Z_1,R_1,ZPt_1,0,10,hArea_1,JetEta_1,D0Eta_1,JetPhi_1,D0Phi_1,D0JetEta_1,Eta_1,JetPt_R_1);
  sumw1+=fillHistos("../sPlot/sWeight_PT_1_10_Cen_0_10.root",D0Pt_1,JetPt_1,Z_1,R_1,ZPt_1,0,10,hArea_1,JetEta_1,D0Eta_1,JetPhi_1,D0Phi_1,D0JetEta_1,Eta_1,JetPt_R_1);

  // sumw2+=fillHistos("../sPlot/sWeight_PT_1_3_Cen_10_20.root",D0Pt_2,JetPt_2,Z_2,R_2,ZPt_2,10,20,hArea_2,JetEta_2,D0Eta_2,JetPhi_2,D0Phi_2,D0JetEta_2,Eta_2,JetPt_R_2);
  // sumw2+=fillHistos("../sPlot/sWeight_PT_3_10_Cen_10_20.root",D0Pt_2,JetPt_2,Z_2,R_2,ZPt_2,10,20,hArea_2,JetEta_2,D0Eta_2,JetPhi_2,D0Phi_2,D0JetEta_2,Eta_2,JetPt_R_2);
  sumw2+=fillHistos("../sPlot/sWeight_PT_1_10_Cen_10_20.root",D0Pt_2,JetPt_2,Z_2,R_2,ZPt_2,10,20,hArea_2,JetEta_2,D0Eta_2,JetPhi_2,D0Phi_2,D0JetEta_2,Eta_2,JetPt_R_2);

  // sumw3+=fillHistos("../sPlot/sWeight_PT_1_3_Cen_20_30.root",D0Pt_3,JetPt_3,Z_3,R_3,ZPt_3,20,30,hArea_3,JetEta_3,D0Eta_3,JetPhi_3,D0Phi_3,D0JetEta_3,Eta_3,JetPt_R_3);
  // sumw3+=fillHistos("../sPlot/sWeight_PT_3_10_Cen_20_30.root",D0Pt_3,JetPt_3,Z_3,R_3,ZPt_3,20,30,hArea_3,JetEta_3,D0Eta_3,JetPhi_3,D0Phi_3,D0JetEta_3,Eta_3,JetPt_R_3);
  sumw3+=fillHistos("../sPlot/sWeight_PT_1_10_Cen_20_30.root",D0Pt_3,JetPt_3,Z_3,R_3,ZPt_3,20,30,hArea_3,JetEta_3,D0Eta_3,JetPhi_3,D0Phi_3,D0JetEta_3,Eta_3,JetPt_R_3);

  // sumw4+=fillHistos("../sPlot/sWeight_PT_1_3_Cen_30_40.root",D0Pt_4,JetPt_4,Z_4,R_4,ZPt_4,30,40,hArea_4,JetEta_4,D0Eta_4,JetPhi_4,D0Phi_4,D0JetEta_4,Eta_4,JetPt_R_4);
  // sumw4+=fillHistos("../sPlot/sWeight_PT_3_10_Cen_30_40.root",D0Pt_4,JetPt_4,Z_4,R_4,ZPt_4,30,40,hArea_4,JetEta_4,D0Eta_4,JetPhi_4,D0Phi_4,D0JetEta_4,Eta_4,JetPt_R_4);
  sumw4+=fillHistos("../sPlot/sWeight_PT_1_10_Cen_30_40.root",D0Pt_4,JetPt_4,Z_4,R_4,ZPt_4,30,40,hArea_4,JetEta_4,D0Eta_4,JetPhi_4,D0Phi_4,D0JetEta_4,Eta_4,JetPt_R_4);

  // sumw5+=fillHistos2("../sPlot/sWeight_PT_1_3_Cen_40_80.root",D0Pt_5,JetPt_5,Z_5,R_5,ZPt_5,D0JPt_5,Hard2D2_5,40,80,hArea_5,JetEta_5,D0Eta_5,JetPhi_5,D0Phi_5,D0JetEta_5,Eta_5,JetPt_R_5);
  // sumw5+=fillHistos2("../sPlot/sWeight_PT_3_10_Cen_40_80.root",D0Pt_5,JetPt_5,Z_5,R_5,ZPt_5,D0JPt_5,Hard2D2_5,40,80,hArea_5,JetEta_5,D0Eta_5,JetPhi_5,D0Phi_5,D0JetEta_5,Eta_5,JetPt_R_5);
  sumw5+=fillHistos2("../sPlot/sWeight_PT_1_10_Cen_40_80.root",D0Pt_5,JetPt_5,Z_5,R_5,ZPt_5,D0JPt_5,Hard2D2_5,40,80,hArea_5,JetEta_5,D0Eta_5,JetPhi_5,D0Phi_5,D0JetEta_5,Eta_5,JetPt_R_5);


  // sumw6+=fillHistos("../sPlot/sWeight_PT_1_3_Cen_10_20.root",D0Pt_6,JetPt_6,Z_6,R_6,ZPt_6,10,40,hArea_6,JetEta_6,D0Eta_6,JetPhi_6,D0Phi_6,D0JetEta_6,Eta_6,JetPt_R_6);
  // sumw6+=fillHistos("../sPlot/sWeight_PT_3_10_Cen_10_20.root",D0Pt_6,JetPt_6,Z_6,R_6,ZPt_6,10,40,hArea_6,JetEta_6,D0Eta_6,JetPhi_6,D0Phi_6,D0JetEta_6,Eta_6,JetPt_R_6);
  sumw6+=fillHistos("../sPlot/sWeight_PT_1_10_Cen_10_20.root",D0Pt_6,JetPt_6,Z_6,R_6,ZPt_6,10,40,hArea_6,JetEta_6,D0Eta_6,JetPhi_6,D0Phi_6,D0JetEta_6,Eta_6,JetPt_R_6);
  // sumw6+=fillHistos("../sPlot/sWeight_PT_1_3_Cen_20_30.root",D0Pt_6,JetPt_6,Z_6,R_6,ZPt_6,10,40,hArea_6,JetEta_6,D0Eta_6,JetPhi_6,D0Phi_6,D0JetEta_6,Eta_6,JetPt_R_6);
  // sumw6+=fillHistos("../sPlot/sWeight_PT_3_10_Cen_20_30.root",D0Pt_6,JetPt_6,Z_6,R_6,ZPt_6,10,40,hArea_6,JetEta_6,D0Eta_6,JetPhi_6,D0Phi_6,D0JetEta_6,Eta_6,JetPt_R_6);
  sumw6+=fillHistos("../sPlot/sWeight_PT_1_10_Cen_20_30.root",D0Pt_6,JetPt_6,Z_6,R_6,ZPt_6,10,40,hArea_6,JetEta_6,D0Eta_6,JetPhi_6,D0Phi_6,D0JetEta_6,Eta_6,JetPt_R_6);
  // sumw6+=fillHistos("../sPlot/sWeight_PT_1_3_Cen_30_40.root",D0Pt_6,JetPt_6,Z_6,R_6,ZPt_6,10,40,hArea_6,JetEta_6,D0Eta_6,JetPhi_6,D0Phi_6,D0JetEta_6,Eta_6,JetPt_R_6);
  // sumw6+=fillHistos("../sPlot/sWeight_PT_3_10_Cen_30_40.root",D0Pt_6,JetPt_6,Z_6,R_6,ZPt_6,10,40,hArea_6,JetEta_6,D0Eta_6,JetPhi_6,D0Phi_6,D0JetEta_6,Eta_6,JetPt_R_6);
  sumw6+=fillHistos("../sPlot/sWeight_PT_1_10_Cen_30_40.root",D0Pt_6,JetPt_6,Z_6,R_6,ZPt_6,10,40,hArea_6,JetEta_6,D0Eta_6,JetPhi_6,D0Phi_6,D0JetEta_6,Eta_6,JetPt_R_6);
  
  // sumw7+=fillHistos("../sPlot/sWeight_PT_1_3_Cen_40_80.root",D0Pt_7,JetPt_7,Z_7,R_7,ZPt_7,40,60,hArea_7,JetEta_7,D0Eta_7,JetPhi_7,D0Phi_7,D0JetEta_7,Eta_7,JetPt_R_7);
  // sumw7+=fillHistos("../sPlot/sWeight_PT_3_10_Cen_40_80.root",D0Pt_7,JetPt_7,Z_7,R_7,ZPt_7,40,60,hArea_7,JetEta_7,D0Eta_7,JetPhi_7,D0Phi_7,D0JetEta_7,Eta_7,JetPt_R_7);
  sumw7+=fillHistos("../sPlot/sWeight_PT_1_10_Cen_40_80.root",D0Pt_7,JetPt_7,Z_7,R_7,ZPt_7,40,60,hArea_7,JetEta_7,D0Eta_7,JetPhi_7,D0Phi_7,D0JetEta_7,Eta_7,JetPt_R_7);

  sumw8+=fillHistos("../sPlot/sWeight_PT_1_10_Cen_0_10.root",D0Pt_8,JetPt_8,Z_8,R_8,ZPt_8,0,10,hArea_8,JetEta_8,D0Eta_8,JetPhi_8,D0Phi_8,D0JetEta_8,Eta_8,JetPt_R_8);

  sumw9+=fillHistos("../sPlot/sWeight_PT_1_10_Cen_10_20.root",D0Pt_9,JetPt_9,Z_9,R_9,ZPt_9,10,40,hArea_9,JetEta_9,D0Eta_9,JetPhi_9,D0Phi_9,D0JetEta_9,Eta_9,JetPt_R_9);
  sumw9+=fillHistos("../sPlot/sWeight_PT_1_10_Cen_20_30.root",D0Pt_9,JetPt_9,Z_9,R_9,ZPt_9,10,40,hArea_9,JetEta_9,D0Eta_9,JetPhi_9,D0Phi_9,D0JetEta_9,Eta_9,JetPt_R_9);
  sumw9+=fillHistos("../sPlot/sWeight_PT_1_10_Cen_30_40.root",D0Pt_9,JetPt_9,Z_9,R_9,ZPt_9,10,40,hArea_9,JetEta_9,D0Eta_9,JetPhi_9,D0Phi_9,D0JetEta_9,Eta_9,JetPt_R_9);

  sumw10+=fillHistos("../sPlot/sWeight_PT_1_10_Cen_40_80.root",D0Pt_10,JetPt_10,Z_10,R_10,ZPt_10,40,80,hArea_10,JetEta_10,D0Eta_10,JetPhi_10,D0Phi_10,D0JetEta_10,Eta_10,JetPt_R_10);

  if(!save){
      DN(D0Pt,14,4,sumw);
      DN(D0Pt_1,2,21,sumw1);
      DN(D0Pt_2,9,22,sumw2);
      DN(D0Pt_3,32,23,sumw3);
      DN(D0Pt_4,6,24,sumw4);   
      DN(D0Pt_5,1,8,sumw5);         

      DN(JetPt,14,4,1);//sumw);
      DN(JetPt_1,2,21,1);//sumw1);
      DN(JetPt_2,9,22,1);//sumw2);
      DN(JetPt_3,32,23,1);//sumw3);
      DN(JetPt_4,6,24,1);//sumw4); 
      DN(JetPt_5,1,8,1);//sumw5);
      
      DN(R,14,4,sumw);
      DN(R_1,2,21,sumw1);
      DN(R_2,9,22,sumw2);
      DN(R_3,32,23,sumw3);
      DN(R_4,6,24,sumw4);
      DN(R_5,1,8,sumw5);
      
      DN(Z,14,4,sumw);
      DN(Z_1,2,21,sumw1);
      DN(Z_2,9,22,sumw2);
      DN(Z_3,32,23,sumw3);
      DN(Z_4,6,24,sumw4);
      DN(Z_5,1,8,sumw5);
  }
  
  D0Pt->GetXaxis()->SetTitle("D^{0} p_{T} (GeV/#it{c})");
  JetPt->GetXaxis()->SetTitle("Jet p_{T} (GeV/#it{c})");
  Z->GetXaxis()->SetTitle("#it{z}");
  R->GetXaxis()->SetTitle("#DeltaR");

  D0Pt->GetYaxis()->SetTitle("dN/d#it{p}_{T} (#it{c}/GeV)");
  JetPt->GetYaxis()->SetTitle("dN/d#it{p}_{T} (#it{c}/GeV)");
  Z->GetYaxis()->SetTitle("(1/N_{jet}) dN/d#it{z}");
  R->GetYaxis()->SetTitle("(1/N_{jet}) ddN/d#DeltaR");

  D0Pt->GetXaxis()->SetRangeUser(0,10);
  JetPt->GetXaxis()->SetRangeUser(3,40);
  R->GetXaxis()->SetRangeUser(0,0.4);
  Z->GetXaxis()->SetRangeUser(0,1);

  D0Pt->GetYaxis()->SetRangeUser(0.00001,D0Pt->GetMaximum()*10.5);
  JetPt->GetYaxis()->SetRangeUser(0.00001,JetPt->GetMaximum()*10.5);
  R->GetYaxis()->SetRangeUser(0.1,R->GetMaximum()*10.5);
  Z->GetYaxis()->SetRangeUser(0.001,Z->GetMaximum()*10.5);
  
  if(!SYS){
      TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
      leg->SetTextSize(0.04);
      leg->AddEntry(R,"0-80%","PE");
      leg->AddEntry(R_1,"0-10%","PE");
      leg->AddEntry(R_2,"10-20%","PE");
      leg->AddEntry(R_3,"20-30%","PE");
      leg->AddEntry(R_4,"30-40%","PE");
      leg->AddEntry(R_5,"40-80%","PE");
      TCanvas *c11 = new TCanvas("c11","c11");
      D0Pt->DrawClone("PE");
      D0Pt_1->DrawClone("PE same");
      D0Pt_2->DrawClone("PE same");
      D0Pt_3->DrawClone("PE same");
      D0Pt_4->DrawClone("PE same");
      D0Pt_5->DrawClone("PE same");
      leg->Draw("same");
      gPad->SetLogy();
      
      TCanvas *c21 =new TCanvas("c21","c21");
      JetPt->DrawClone("PE");
      // TCanvas *c211 =new TCanvas("c211","c211");
      JetPt_1->DrawClone("PE same");
      //TCanvas *c212 =new TCanvas("c212","c212");
      JetPt_2->DrawClone("PE same");
      //TCanvas *c213 =new TCanvas("c213","c213");
      JetPt_3->DrawClone("PE same");
      //TCanvas *c21 =new TCanvas("c24","c24");
      JetPt_4->DrawClone("PE same");
      JetPt_5->DrawClone("PE same");
      leg->Draw("same");
      gPad->SetLogy();
      TCanvas *c31 =new TCanvas("c31","c31");
      Z->DrawClone("PE");
      Z_1->DrawClone("PE same");
      Z_2->DrawClone("PE same");
      Z_3->DrawClone("PE same");
      Z_4->DrawClone("PE same");
      Z_5->DrawClone("PE same");
      leg->Draw("same");
      gPad->SetLogy();
      TCanvas *c41 =new TCanvas("c41","c41");
      R->DrawClone("PE");
      R_1->DrawClone("PE same");
      R_2->DrawClone("PE same");
      R_3->DrawClone("PE same");
      R_4->DrawClone("PE same");
      R_5->DrawClone("PE same");
      leg->Draw("same");
      gPad->SetLogy();
      
      PT2D->GetXaxis()->SetTitle("Jet p^{corr.}_{T} (GeV/c)");
      PT2D->GetYaxis()->SetTitle("D^{0} p_{T} (GeV/c)");
      PT2D->GetZaxis()->SetRangeUser(2,PT2D->GetMaximum()*5);
      TCanvas *c411 =new TCanvas("c411","c411");
      PT2D->Draw("COLZ");
  }
  if(save){
      TFile *outf = new TFile(NAME,"RECREATE");
      centrality->Write();
      D0Pt->Write("D0Pt");
      cout << "Integral here " << D0Pt_1->Integral() << " " << sumw1 << endl;
      D0Pt_8->Write("D0Pt_0_10");
      D0Pt_2->Write("D0Pt_10_20");
      D0Pt_3->Write("D0Pt_20_30");
      D0Pt_4->Write("D0Pt_30_40"); 
      D0Pt_10->Write("D0Pt_40_80");
      D0Pt_9->Write("D0Pt_10_40");
      D0JPt_5->Write("D0JPt_40_80");
      JetPt->Write("JetPt");
      JetPt_8->Write("JetPt_0_10");
      JetPt_2->Write("JetPt_10_20");
      JetPt_3->Write("JetPt_20_30");
      JetPt_4->Write("JetPt_30_40");
      JetPt_10->Write("JetPt_40_80");
      JetPt_9->Write("JetPt_10_40");
      Z->Write("Z");
      Z_8->Write("Z_0_10");
      Z_2->Write("Z_10_20");
      Z_3->Write("Z_20_30");
      Z_4->Write("Z_30_40");
      Z_10->Write("Z_40_80");
      Z_9->Write("Z_10_40");
      ZPt->Write("ZPt");
      ZPt_8->Write("ZPt_0_10");
      ZPt_2->Write("ZPt_10_20");
      ZPt_3->Write("ZPt_20_30");
      ZPt_4->Write("ZPt_30_40");
      ZPt_10->Write("ZPt_40_80");
      ZPt_9->Write("ZPt_10_40");
      R->Write("R");
      R_8->Write("R_0_10");
      R_2->Write("R_10_20");
      R_3->Write("R_20_30");
      R_4->Write("R_30_40");
      R_10->Write("R_40_80");
      R_9->Write("R_10_40");
      
      D0Pt_7->Write("D0Pt_40_60");
      JetPt_7->Write("JetPt_40_60");
      Z_7->Write("Z_40_60");
      ZPt_7->Write("ZPt_40_60");
      R_7->Write("R_40_60");
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
      

      JetPt_R->Write("JetPt_R");
      JetPt_R_1->Write("JetPt_R_0_10");
      JetPt_R_2->Write("JetPt_R_10_20");
      JetPt_R_3->Write("JetPt_R_20_30");
      JetPt_R_4->Write("JetPt_R_30_40");
      JetPt_R_5->Write("JetPt_R_40_80");
      JetPt_R_6->Write("JetPt_R_10_40");
      JetPt_R_7->Write("JetPt_R_40_60");
  }

}

void doSys(){
    SYS = true;
    fillSys(D0Sys_010,1);
    fillSys(D0Sys_1040,2);
    fillSys(D0Sys_4080,3);
    for(int i = 0;i<166;i++){
      cout << ">>> On systematic loop " << i+1 << endl;
      sprintf(NAME,"systematics_new_ratios/Histograms3_D05GeV_%i.root",i);
      cout << " >> Will save histograms as " << NAME << endl;
      shuffleEff();
      applyWeightsNew(1);
    }

}
