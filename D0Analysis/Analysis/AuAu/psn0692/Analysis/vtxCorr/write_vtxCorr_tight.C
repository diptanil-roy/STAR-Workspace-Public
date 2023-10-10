#include "../anaCuts.h"
#include "myFunction.h"
void write_vtxCorr_tight() {
    char dir[250];
    char CMD[250];
    char name[250];
    char title[250];
    
    sprintf(dir,"data");
    sprintf(CMD,"[ -d %s ] || mkdir -p %s",dir,dir);
    gSystem->Exec(CMD);
    
    const int nType = 2;
    const char nameDir[nType][100] = {"./out", "./out"};
    const char nameType[nType][100] = {"", "_sys"};
    const char nameType1[nType][100] = {"gRef cut 1", "gRef cut 2"};
    const char nameType2[nType][100] = {"tight", "sys"};
    
    //read eff from pure hijing
    TH2F* h2Base_hj[nType];
    TH2F* h2Eff_hj[nType];
    for(int i=0; i<nType; i++) {
        TFile* fin1 = new TFile(Form("%s/effPureHijing%s.root",nameDir,nameType[i]));
        h2Base_hj[i] = (TH2F*)fin1->Get("h2McPtCent");
        h2Base_hj[i]->SetDirectory(0);
        h2Eff_hj[i] = (TH2F*)fin1->Get("h2PtCent_tight_pt1");
        h2Eff_hj[i]->SetDirectory(0);
        fin1->Close();
    }
    
    //read eff from fast simu based on hijing
    TH2F* h2Base_fs[nType];
    TH2F* h2Eff_fs[nType];
    for(int i=0; i<nType; i++) {
        TFile* fin2 = new TFile(Form("%s/eff.FastSimu.Hijing%s.root",nameDir,nameType[i]));
        h2Base_fs[i] = (TH2F*)fin2->Get("h2McPtCent");
        h2Base_fs[i]->SetDirectory(0);
        h2Eff_fs[i] = (TH2F*)fin2->Get("h2PtCent_tight_pt1");
        h2Eff_fs[i]->SetDirectory(0);
        fin2->Close();
    }
    
    //calculate eff
    TH1F* heffHj[nType][ncent_vtx];
    TH1F* heffFs[nType][ncent_vtx];
    for(int icent=0; icent<ncent_vtx; icent++) {
        for(int itype=0; itype<nType; itype++) {
            TH1F* hbase_hj = (TH1F*)h2Base_hj[itype]->ProjectionX(Form("hbaseHj_%i_%i",icent,itype),icent+1,icent+1)->Rebin(npt,Form("hbaseHjRebin_%i_%i",icent,itype),nptbin);
            heffHj[itype][icent] = (TH1F*)h2Eff_hj[itype]->ProjectionX(Form("heffHj_%i",icent),icent+1,icent+1)->Rebin(npt,Form("heffHjRebin_%i_%i",icent,itype),nptbin);
            //heffHj[itype][icent]->Divide(hbase_hj);
            heffHj[itype][icent] = histDivide(heffHj[itype][icent],hbase_hj,true);
            
            TH1F* hbase_fs = (TH1F*)h2Base_fs[itype]->ProjectionX(Form("hbaseFs_%i_%i",icent,itype),icent+1,icent+1)->Rebin(npt,Form("hbaseFsRebin_%i_%i",icent,itype),nptbin);
            heffFs[itype][icent] = (TH1F*)h2Eff_fs[itype]->ProjectionX(Form("heffFs_%i_%i",icent,itype),icent+1,icent+1)->Rebin(npt,Form("heffFsRebin_%i_%i",icent,itype),nptbin);
            //heffFs[itype][icent]->Divide(hbase_fs);
            heffFs[itype][icent] = histDivide(heffFs[itype][icent],hbase_fs,true);
        }
    }
    
    // cal. ratio (fast sim./ pure hij.)
    TH1F* hRatio[nType][ncent_vtx];
    for(int icent=0; icent<ncent_vtx; icent++) {
        for(int itype=0; itype<nType; itype++) {
            hRatio[itype][icent] = (TH1F*)heffHj[itype][icent]->Clone(Form("hRatio_%i_%i",itype,icent));
            hRatio[itype][icent]->Divide(heffFs[itype][icent]);
        }
    }
    
    // all cent bin using data by data
    // cal. correction ratio and its sys error: 0-80%
    /*TH1F* hmean[ncentAll_vtx];
     TH1F* herr[ncentAll_vtx];
     for(int icent=0; icent<ncent_vtx; icent++) {
     cout << Form("mean_%s",nameCentAll_vtx1[icent]) << endl;
     hmean[icent] = (TH1F*)hRatio[0][icent]->Clone(Form("mean_%s",nameCentAll_vtx1[icent]));
     herr[icent] = (TH1F*)hRatio[0][icent]->Clone(Form("sysError_%s",nameCentAll_vtx1[icent]));
     for(int ibin=1; ibin<=hRatio[0][icent]->GetNbinsX(); ibin++) {
     float num1 = hRatio[0][icent]->GetBinContent(ibin);
     float num2 = hRatio[1][icent]->GetBinContent(ibin);
     float error;
     if((num1+num2)==0) error = 0;
     else error = fabs(num1 - num2)/(num1 + num2);
     herr[icent]->SetBinContent(ibin,error);
     hmean[icent]->SetBinContent(ibin,(num1+num2)/2.);
     herr[icent]->SetBinError(ibin,0);
     hmean[icent]->SetBinError(ibin,fabs(num1 - num2)/2.);
     }
     }
     
     // cal. correction ratio and its sys error: 0-10%
     for(int icent=ncent_vtx; icent<ncentAll_vtx; icent++) {
     cout << Form("mean_%s",nameCentAll_vtx1[icent]) << endl;
     hmean[icent] = (TH1F*)hmean[ncent_vtx-1]->Clone(Form("mean_%s",nameCentAll_vtx1[icent]));
     herr[icent] = (TH1F*)herr[ncent_vtx-1]->Clone(Form("sysError_%s",nameCentAll_vtx1[icent]));
     }*/
    
    // 20-80% using data by data, 0-20% meas value set to be 1
    // cal. correction ratio and its sys error: 20-80%
    TH1F* hmean[ncentAll_vtx];
    TH1F* herr[ncentAll_vtx];
    int centCut = -1;
    if(ncent_vtx>6) centCut = 6;
    else centCut = ncent_vtx;
    for(int icent=0; icent<centCut; icent++) {
        cout << Form("mean_%s",nameCentAll_vtx1[icent]) << endl;
        hmean[icent] = (TH1F*)hRatio[0][icent]->Clone(Form("mean_%s",nameCentAll_vtx1[icent]));
        herr[icent] = (TH1F*)hRatio[0][icent]->Clone(Form("sysError_%s",nameCentAll_vtx1[icent]));
        for(int ibin=1; ibin<=hRatio[0][icent]->GetNbinsX(); ibin++) {
            float num1 = hRatio[0][icent]->GetBinContent(ibin);
            float num2 = hRatio[1][icent]->GetBinContent(ibin);
            float error;
            if((num1+num2)==0) error = 0;
            else error = fabs(num1 - num2)/(num1 + num2);
            error = sqrt(pow(error,2) + pow(0.05,2)); // 5% additional systematic uncertainties was assigned
            herr[icent]->SetBinContent(ibin,error);
            // hmean[icent]->SetBinContent(ibin,(num1+num2)/2.);
            hmean[icent]->SetBinContent(ibin,num2);
            herr[icent]->SetBinError(ibin,0);
            // hmean[icent]->SetBinError(ibin,fabs(num1 - num2)/2.);
            hmean[icent]->SetBinError(ibin,num2*error);
        }
    }
    
    // cal. correction ratio and its sys error: 0-10%
    for(int icent=centCut; icent<ncentAll_vtx-2; icent++) {
        cout << Form("mean_%s",nameCentAll_vtx1[icent]) << endl;
        hmean[icent] = (TH1F*)hRatio[0][icent]->Clone(Form("mean_%s",nameCentAll_vtx1[icent]));
        herr[icent] = (TH1F*)hRatio[0][icent]->Clone(Form("sysError_%s",nameCentAll_vtx1[icent]));
        for(int ibin=1; ibin<=hmean[icent]->GetNbinsX(); ibin++) {
            //float num1 = 1;
            //float num2 = 1;
            //float error;
            //if((num1+num2)==0) error = 0;
            //else error = fabs(num1 - num2)/(num1 + num2);
            //herr[icent]->SetBinContent(ibin,error);
            //hmean[icent]->SetBinContent(ibin,(num1+num2)/2.);
            float num = 1;
            float error = 0.05;
            // float num1 = hRatio[0][icent]->GetBinContent(ibin);
            // float num2 = hRatio[1][icent]->GetBinContent(ibin);
            // float error;
            // if((num1+num2)==0) error = 0;
            // else error = fabs(num1 - num2)/(num1 + num2);
            // error = sqrt(pow(error,2) + pow(0.05,2)); // 5% additional systematic uncertainties was assigned
            herr[icent]->SetBinContent(ibin,error);
            hmean[icent]->SetBinContent(ibin,num);
            // hmean[icent]->SetBinContent(ibin,num2);
            herr[icent]->SetBinError(ibin,0);
            hmean[icent]->SetBinError(ibin,error*num);
            // hmean[icent]->SetBinError(ibin,error*num2);
        }
    }
    
    // cal. correction ratio and its sys error: 0-20%
    for(int icent=centCut+1; icent<ncentAll_vtx; icent++) {
        cout << Form("mean_%s",nameCentAll_vtx1[icent]) << endl;
        hmean[icent] = (TH1F*)hRatio[0][ncent_vtx-1]->Clone(Form("mean_%s",nameCentAll_vtx1[icent]));
        herr[icent] = (TH1F*)hRatio[0][ncent_vtx-1]->Clone(Form("sysError_%s",nameCentAll_vtx1[icent]));
        for(int ibin=1; ibin<=hmean[icent]->GetNbinsX(); ibin++) {
            //float num1 = 1;
            //float num2 = 1;
            //float error;
            //if((num1+num2)==0) error = 0;
            //else error = fabs(num1 - num2)/(num1 + num2);
            //herr[icent]->SetBinContent(ibin,error);
            //hmean[icent]->SetBinContent(ibin,(num1+num2)/2.);
            // float num = 1;
            // float error = 0.05;
            float num1 = hRatio[0][ncent_vtx-1]->GetBinContent(ibin);
            float num2 = hRatio[1][ncent_vtx-1]->GetBinContent(ibin);
            float error;
            if((num1+num2)==0) error = 0;
            else error = fabs(num1 - num2)/(num1 + num2);
            error = sqrt(pow(error,2) + pow(0.05,2)); // 5% additional systematic uncertainties was assigned
            herr[icent]->SetBinContent(ibin,error);
            // hmean[icent]->SetBinContent(ibin,num);
            hmean[icent]->SetBinContent(ibin,num2);
            herr[icent]->SetBinError(ibin,0);
            // hmean[icent]->SetBinError(ibin,error*num);
            hmean[icent]->SetBinError(ibin,error*num2);
        }
    }
    
    // write
    sprintf(name,"%s/vtxCorr_%s.root",dir,nameType2[0]);
    TFile* fout = new TFile(name,"RECREATE");
    for(int icent=0; icent<ncentAll_vtx; icent++) {
        hmean[icent]->Write();
        herr[icent]->Write();
    }
    fout->Close();
}
