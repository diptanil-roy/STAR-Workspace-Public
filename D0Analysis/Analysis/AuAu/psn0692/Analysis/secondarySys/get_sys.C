#include "../anaCuts.h"
void get_sys() {
    ifstream in;
    ofstream out;
    float y[npt], yerr[npt], ybase[npt], yerrbase[npt];
    float sys;
    
    TFile* fin = new TFile("../secondarySys/data/eff.root","open");
    TH1D* heffBinD0[ncent];
    for(int icent=0; icent<ncent; icent++) {
        heffBinD0[icent] = (TH1D*)fin->Get(Form("heffBinD0_%i",icent));
        heffBinD0[icent]->SetDirectory(0);
    }
    fin->Close();

    TFile* finbase = new TFile("../default/data/eff.root","open");
    TH1D* heffBinD0base[ncent];
    for(int icent=0; icent<ncent; icent++) {
        heffBinD0base[icent] = (TH1D*)finbase->Get(Form("heffBinD0_%i",icent));
        heffBinD0base[icent]->SetDirectory(0);
    }
    finbase->Close();

    for(int icent=0; icent<ncent; icent++) {
        out.open(Form("data/yieldSys_%s.txt",nameCent1[icent]));
        for(int ipt=0; ipt<npt; ipt++) {
            float pt = 0.5*(nptbin[ipt]+nptbin[ipt+1]);
            int bin = heffBinD0[icent]->FindBin(pt);
            float eff = heffBinD0[icent]->GetBinContent(bin);
            float effbase = heffBinD0base[icent]->GetBinContent(bin);
            sys = fabs(eff-effbase)/effbase;
            out << ipt << "\t" << sys << endl;
        }
        out.close();
    }

}
