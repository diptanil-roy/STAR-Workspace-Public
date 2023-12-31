#include "../anaCuts.h"
#include "myFunction.h"
#include "../myConst.h"
void plot_statErr() {
    globalSetting();
    char dir[250];
    char name[250];
    char CMD[250];
    char lname[16][250];
    //TPaveStats* ptstates;
    TLegend* legend;
    TH1F* h0;
    
    sprintf(dir,"pic/stat");
    sprintf(CMD,"[ -d %s ] || mkdir -p %s",dir,dir);
    gSystem->Exec(CMD);
    
    const int nType = 3;
    const char nameDir[nType][100] = {"./out", "./out", "./out"};
    const char nameType[nType][100] = {"", "_half1", "_half2"};
    const char nameType1[nType][100] = {"default", "half sample1", "half sample2"};
    
    //read eff from pure hijing
    TH2F* h2Base_hj[nType];
    TH2F* h2Eff_hj[nType];
    for(int i=0; i<nType; i++) {
        TFile* fin1 = new TFile(Form("%s/effPureHijing%s.root",nameDir,nameType[i]));
        h2Base_hj[i] = (TH2F*)fin1->Get("h2McPtCent");
        h2Base_hj[i]->SetDirectory(0);
        h2Eff_hj[i] = (TH2F*)fin1->Get("h2PtCent_pt1");
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
        h2Eff_fs[i] = (TH2F*)fin2->Get("h2PtCent_pt1");
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
    
    
    // cal. error
    TH1F* herr[ncent_vtx];
    for(int icent=0; icent<ncent_vtx; icent++) {
        herr[icent] = (TH1F*)hRatio[0][icent]->Clone(Form("statError_%s",nameCent_vtx1[icent]));
        for(int ibin=1; ibin<=hRatio[0][icent]->GetNbinsX(); ibin++) {
            float base = hRatio[0][icent]->GetBinContent(ibin);
            float num1 = hRatio[1][icent]->GetBinContent(ibin);
            float num2 = hRatio[2][icent]->GetBinContent(ibin);
            float error = (fabs(num1-base) + fabs(num2-base))/2.0;
            //if(fabs(num2-base)>error) error = fabs(num2-base);
            if(base==0) error = 0;
            else error /= base;
            herr[icent]->SetBinContent(ibin,error);
        }
    }
    
    // cal. error directly from statistics
    TH1F* herr2[ncent_vtx];
    for(int icent=0; icent<ncent_vtx; icent++) {
        herr2[icent] = (TH1F*)hRatio[0][icent]->Clone(Form("statError2_%s",nameCent_vtx1[icent]));
        for(int ibin=1; ibin<=hRatio[0][icent]->GetNbinsX(); ibin++) {
            float base = hRatio[0][icent]->GetBinContent(ibin);
            // float num1 = hRatio[1][icent]->GetBinContent(ibin);
            // float num2 = hRatio[2][icent]->GetBinContent(ibin);
            // float error = (fabs(num1-base) + fabs(num2-base))/2.0;
            float error = hRatio[0][icent]->GetBinError(ibin);
            //if(fabs(num2-base)>error) error = fabs(num2-base);
            if(base==0) error = 0;
            else error /= base;
            herr2[icent]->SetBinContent(ibin,error);
        }
    }
    
    // plot
    TCanvas* c1 = new TCanvas("c1", "A Canvas",10,10,800,800);
    setPad(c1);
    const float linewidth = 2;
    const float markersize =  1.5;
    for(int icent=0; icent<ncent_vtx; icent++) {
        for(int itype=0; itype<nType; itype++) {
            hRatio[itype][icent]->SetMarkerStyle(MARKERSTYLE[itype]);
            hRatio[itype][icent]->SetMarkerColor(COLOR[itype]);
            hRatio[itype][icent]->SetMarkerSize(markersize);
            hRatio[itype][icent]->SetLineColor(COLOR[itype]);
            hRatio[itype][icent]->SetLineWidth(linewidth);
        }
        herr[icent]->SetLineWidth(linewidth);
        herr[icent]->SetLineColor(COLOR[nType]);
        herr[icent]->SetMarkerStyle(24);
        herr[icent]->SetMarkerColor(COLOR[nType]);
        herr2[icent]->SetLineWidth(linewidth);
        herr2[icent]->SetLineColor(COLOR[nType-1]);
        herr2[icent]->SetMarkerStyle(24);
        herr2[icent]->SetMarkerColor(COLOR[nType-1]);
        
        for(int itype=0; itype<nType; itype++) {
            if(itype==0) {
                hRatio[itype][icent]->Draw("e");
                setHisto(hRatio[itype][icent],"","p_{T} (GeV/c)", "Eff. Ratio (pure HIJING / fast Sim)");
                hRatio[itype][icent]->GetYaxis()->SetRangeUser(0,1.8);
                //hRatio[itype][icent]->GetYaxis()->SetRangeUser(1e-3,9);
            }
            else hRatio[itype][icent]->Draw("esame");
        }
        
        herr[icent]->Draw("histsame");
        herr2[icent]->Draw("histsame");
        
        if(!legend) {
            const float legy_up = 0.88;
            const float legy_lw = legy_up-(nType+1)*0.05;
            legend = new TLegend(0.2,legy_lw,0.6,legy_up);
            legend->SetFillStyle(0);
            legend->SetFillColor(10);
            legend->SetBorderSize(0);
            legend->SetTextSize(0.05);
            legend->SetTextFont(132);
            for(int itype=0; itype<nType; itype++) legend->AddEntry(hRatio[itype][icent],nameType1[itype],"p");
            legend->AddEntry(herr[icent],"error half-stat","l");
            legend->AddEntry(herr2[icent],"error full sample","l");
        }
        legend->Draw();
        
        sprintf(name,"AuAu @200GeV %s", nameCent_vtx[icent]);
        drawLatex(0.2,0.89,name,132,0.05,1);
        
        //gPad->SetLogy();
        sprintf(name,"%s/halfStat_%s.gif",dir,nameCent_vtx1[icent]);
        c1->SaveAs(name);
    }
}
