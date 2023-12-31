#include "../anaCuts.h"
#include "../myConst.h"
#include "../myFunction.h"
void doPtShift() {
    globalSetting();
    char name[250];
    char CMD[250];
    char dir[250];
    TLegend* legend;
    TH1F* h0;
    TCanvas* c1 = new TCanvas("c1", "A Canvas",10,10,800,800);
    setPad(c1);
    
    sprintf(dir,"pic");
    sprintf(CMD,"[ -d %s ] || mkdir -p %s",dir,dir);
    gSystem->Exec(CMD);
    
    //read spectra without pt shift
    TFile* fin = new TFile("../sys/D0_Spectra_Run14_HFT_beforePtShift.root");
    TGraphErrors* gD0_err[ncent];
    TGraphErrors* gD0_sys[ncent];
    for(int icent=0; icent<ncent; icent++) {
        gD0_err[icent] = (TGraphErrors*)fin->Get(Form("gD0_err_%s",nameCent1[icent]));
        gD0_sys[icent] = (TGraphErrors*)fin->Get(Form("gD0_sys_%s",nameCent1[icent]));
    }
    fin->Close();
    
    //define fit function, now use levy function
    char funcString[200];
    char funcString_time_pt[200];
    double m0 = 1.8645;//D0-1.8645, D+/- - 1.8693;
    sprintf(funcString,"(1/(2*TMath::Pi()))*[0]*([2]-1)*([2]-2)/([2]*[1]*([2]*[1]+%f*([2]-2)))*TMath::Power(([2]*[1]+TMath::Sqrt(x[0]*x[0]+%f*%f)-%f)/([2]*[1]),-[2])",m0,m0,m0,m0); // dN/pTdpTdy
    sprintf(funcString_time_pt,"(1/(2*TMath::Pi()))*[0]*([2]-1)*([2]-2)/([2]*[1]*([2]*[1]+%f*([2]-2)))*TMath::Power(([2]*[1]+TMath::Sqrt(x[0]*x[0]+%f*%f)-%f)/([2]*[1]),-[2])*x[0]",m0,m0,m0,m0); // dN/dpTdy
    TF1 *flevy = new TF1("flevy",funcString,0,20);
    TF1 *flevytmp = new TF1("flevytmp",funcString,0,20);
    // flevy->SetParameters(6.25342e-01,3.34550e-01,1.86042e+01);
    flevy->SetParameters(5.79693e-01,3.47947e-01,1.46732e+01);
    flevy->SetLineColor(kGreen+2);
    // flevy->SetLineStyle(9);
    TF1 *flevy_time_pt = new TF1("flevy_time_pt",funcString_time_pt,0,20);
    
    // Plots 1 : Draw the spectra without any shift
    // ------------------------------------------------
    iCent = 2;// us 20-40% do the shift
    h0 = new TH1F("","",10,0,10);
    setHisto(h0, "", "p_{T} (GeV/c)", "d^{2}N/(N_{ev}2#pip_{T}dp_{T}dy)");
    h0->GetYaxis()->SetRangeUser(0.3*gD0_err[iCent]->GetY()[gD0_err[iCent]->GetN()-1],10*gD0_err[iCent]->GetY()[iCent]);
    h0->GetXaxis()->SetRangeUser(0,10);
    h0->Draw();
    gD0_err[iCent]->SetLineColor(2);
    gD0_err[iCent]->SetMarkerSize(2);
    gD0_err[iCent]->SetMarkerStyle(20);
    gD0_err[iCent]->SetMarkerColor(2);
    gD0_err[iCent]->DrawClone("psame");
    gD0_err[iCent]->Fit(flevy,"INOR","",0,10);
    // gD0_err[iCent]->Fit(flevy,"INOR","",0,10);
    flevy->DrawClone("same");
    gPad->SetLogy();
    sprintf(name,"%s/D0_%s_beforeShift.pdf",dir,nameCent1[iCent]);
    c1->SaveAs(name);
    gPad->SetLogy(1);
    
   // return;

    // plots 2 : Draw all shifts
    // --------------------------
    int NPoint = gD0_err[iCent]->GetN();
    const int nloops = 10;
    Int_t styleindex[6] = {24,25,26,29,27,28};
    Int_t colorindex[10] = {kRed-2,kOrange-2,kSpring-2,kGreen-2,kTeal-2,kCyan-2,kBlue-2,kViolet-2,kMagenta-2,kBlack};
    //h0->Draw();
    ofstream out("ptShift.txt");
    out << "0: ";
    for(int j=0;j<NPoint;j++) out << gD0_err[iCent]->GetX()[j] << "\t";
    out << endl;
    for(int n=0;n<nloops;n++) {
        // flevy->Print();
        flevy_time_pt->SetParameters(flevy->GetParameters());
        // flevy_time_pt->Print();
        
        for(int j=0;j<NPoint;j++)
        {
            Float_t x = gD0_err[iCent]->GetX()[j];
            Float_t y = gD0_err[iCent]->GetY()[j];
            Float_t ey = gD0_err[iCent]->GetEY()[j];
            Float_t x_corr = flevy_time_pt->Integral(nptbin[j],nptbin[j+1])/flevy->Integral(nptbin[j],nptbin[j+1]);
            gD0_err[iCent]->SetPoint(j,x_corr,y);
            gD0_err[iCent]->SetPointError(j,0,ey);
        }
        flevy->SetParameters(6.25342e-01,3.34550e-01,1.86042e+01);
        gD0_err[iCent]->Fit(flevy,"INOR","",0,10);
        gD0_err[iCent]->Fit(flevy,"INOR","",0,10);
        gD0_err[iCent]->SetMarkerColor(colorindex[n]);
        gD0_err[iCent]->SetFillColor(colorindex[n]);
        gD0_err[iCent]->SetLineColor(colorindex[n]);
        flevy->SetLineColor(colorindex[n]);
        gD0_err[iCent]->DrawClone("psame");
        flevy->DrawClone("same");
        out << n+1 << ": ";
        for(int j=0;j<NPoint;j++) out << gD0_err[iCent]->GetX()[j] << "\t";
        out << endl;
    }
    out.close();
    gPad->SetLogy();
    sprintf(name,"%s/D0_%s_ShiftProcess.pdf",dir,nameCent1[iCent]);
    c1->SaveAs(name);
    gPad->SetLogy(1);
    
    // ==== Create a TGraphErrors to hold the confidence intervals
    const Int_t NCL = 200;
    TGraphErrors *grint;
    grint = new TGraphErrors(NCL);
    grint->SetTitle("Fitted line with .68 conf. band");
    for (int i=0; i<NCL; i++) grint->SetPoint(i, 0 + (float)i/20., 0);
    
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint,0.68);
    grint->SetName(Form("flevy_err_band_%s",nameCent1[iCent]));
    
    // ==== Write
    TFile* fout = new TFile("D0_Spectra_Run14HFT.root","RECREATE");
    gD0_err[iCent]->SetMarkerColor(1);
    for(int i=0; i<gD0_err[iCent]->GetN(); i++) gD0_sys[iCent]->GetX()[i] = gD0_err[iCent]->GetX()[i];
    fout->cd();
    gD0_err[iCent]->Write(Form("gD0_err_%s",nameCent1[iCent]));
    gD0_sys[iCent]->Write(Form("gD0_sys_%s",nameCent1[iCent]));
    flevy->Write(Form("flevy_%s",nameCent1[iCent]));
    flevy_time_pt->SetParameters(flevy->GetParameters());
    flevy_time_pt->Write(Form("flevy_time_pt_%s",nameCent1[iCent]));
    grint->Write();

    double tmp[3];
    for(int i=0;i<3;i++) tmp[i] = flevy->GetParameter(i);
    
    // ==== Fit and Write other centrality
    for(icent=0; icent<ncent; icent++) {
      if (icent == iCent) continue;
        Int_t NPoint = gD0_err[iCent]->GetN();
        for(int ip=0;ip<NPoint;ip++)
        {
            gD0_err[icent]->GetX()[ip] = gD0_err[iCent]->GetX()[ip];
            gD0_sys[icent]->GetX()[ip] = gD0_err[iCent]->GetX()[ip];
        }
        
        // flevy->SetParameters(0.00346574, 8.37739, 0.260452);
        if(icent == 0) 
        {
          // flevy->SetParLimits(1, tmp[1]*0.5, tmp[1]*1.5);
          // flevy->SetParLimits(2, tmp[2]*0.5, tmp[2]*1.5);
          flevytmp->SetParameters(6.25342e-01,3.34550e-01,1.86042e+01);
          flevytmp->FixParameter(1, tmp[1]);
          flevytmp->FixParameter(2, tmp[2]);
          for(int ifit=0;ifit<3;ifit++) gD0_err[icent]->Fit(flevytmp,"INOR","",0,10);
          flevy_time_pt->SetParameters(flevytmp->GetParameters());
        }
        else
        {
          flevy->SetParameters(6.25342e-01,3.34550e-01,1.86042e+01);
          for(int ifit=0;ifit<3;ifit++) gD0_err[icent]->Fit(flevy,"INOR","",0,10);
          flevy_time_pt->SetParameters(flevy->GetParameters());
        }


        // for(int ifit=0;ifit<5;ifit++) gD0_err[icent]->Fit(flevy,"INOR","",0,10);
        grint = new TGraphErrors(NCL);
        grint->SetTitle("Fitted line with .68 conf. band");
        for (int i=0; i<NCL; i++) grint->SetPoint(i, 0 + (float)i/20., 0);
        (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint,0.68);
        grint->SetName(Form("flevy_err_band_%s",nameCent1[icent]));
        
        fout->cd();
        gD0_err[icent]->Write(Form("gD0_err_%s",nameCent1[icent]));
        gD0_sys[icent]->Write(Form("gD0_sys_%s",nameCent1[icent]));
        if(icent==0)
        flevytmp->Write(Form("flevy_%s",nameCent1[icent]));
        else 
        flevy->Write(Form("flevy_%s",nameCent1[icent]));
        flevy_time_pt->Write(Form("flevy_time_pt_%s",nameCent1[icent]));
        grint->Write();
    }
    
    fout->Close();
}
