#include "../myFunction.h"
#include "../myConst.h"
void extract_txt() {
    globalSetting();
    char dir[250];
    char name[250];
    char title[250];
    //TString CMD = 0;
    char CMD[250];
    TLegend* legend1;
    TLegend* legend2;
    TH1F* h0;
    
    sprintf(dir,"pic");
    sprintf(CMD,"[ -d %s ] || mkdir -p %s",dir,dir);
    gSystem->Exec(CMD);
    
    const int ncent = 6;
    const char nameCent[ncent][250] = {"0-80%","0-10%", "10-20%", "20-40%", "40-60%", "60-80%"};
    const char nameCent1[ncent][250] = {"0-80% [#times8]","0-10% [#times1]", "10-20% [/3]", "20-40% [/6]", "40-60% [/12]", "60-80% [/24]"};
    const char nameCentXL[ncent][250] = {"0_80","0_10", "10_20", "20_40", "40_60", "60_80"};
    const float scale[ncent] = {8., 1., 1/3., 1/6., 1./12, 1./24};
    
    //Read spectra
    //1. from xiaolong
    TGraphErrors* gD0err_xl[ncent];
    TGraphErrors* gD0sys_xl[ncent];
    TF1* fLevy[ncent];
    TFile* fin1 = new TFile("D0_Spectra_Run14HFT.root");
    for(int icent=0; icent<ncent; icent++) {
        gD0err_xl[icent] = (TGraphErrors*)fin1->Get(Form("gD0_err_%s",nameCentXL[icent]));
        gD0sys_xl[icent] = (TGraphErrors*)fin1->Get(Form("gD0_sys_%s",nameCentXL[icent]));
    }
    fin1->Close();
    
    //extract
    for(int icent=0; icent<ncent; icent++) {
        ofstream out(Form("AuAu_D0_%s.txt",nameCentXL[icent]));
        for(int i=0; i<gD0sys_xl[icent]->GetN(); i++) {
            out << gD0sys_xl[icent]->GetX()[i] << "\t" << gD0sys_xl[icent]->GetY()[i] << "\t" << sqrt(pow(gD0err_xl[icent]->GetEY()[i],2)+pow(gD0sys_xl[icent]->GetEY()[i],2)) << endl;
        }
        out.close();
    }
}
