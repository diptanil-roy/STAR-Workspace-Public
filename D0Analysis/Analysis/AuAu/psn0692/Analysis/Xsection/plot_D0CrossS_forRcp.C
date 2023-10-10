#include "../myFunction.h"
#include "../anaCuts.h"
void plot_D0CrossS_forRcp() {
    globalSetting();
    char dir[250];
    char name[250];
    char title[250];
    //TString CMD = 0;
    char CMD[250];
    TLegend* legend;
    TH1F* h0;
    
    // const int nCent_Draw = 5;
    // const char nameCent_Draw[nCent_Draw][100] = {"60-80%", "40-60%", "20-40%", "10-20%", "0-10%"};
    // const char nameCent1_Draw[nCent_Draw][100] = {"60_80", "40_60", "20_40", "10_20", "0_10"};
    // const float centLw_Draw[ncent] = {0, 2, 4, 6, 7}; // >=    0-80%: 1-9
    // const float centUp_Draw[ncent] = {2, 4, 6, 7, 9}; // <
    // const float Nbin_Draw[nCent_Draw] = {21.37396, 91.37100, 288.35051, 579.89409, 938.80170};
    // const float NbinErr_Draw[nCent_Draw] = {8.93878, 21.05409, 30.39279, 28.86320, 26.28048};
    
    const int nCent_Draw = 10;
    const char nameCent_Draw[nCent_Draw][100] = {"0-10%","10-20%","20-40%","40-60%","60-80%","10-40%","40-80%","0-80%","10-80%","30-50%"};
    const char nameCent1_Draw[nCent_Draw][100] = {"0_10","10_20","20_40","40_60","60_80","10_40","40_80","0_80","10_80","30_50"};
    const float centLw_Draw[ncent] = {7,6,4,2,0,4,0,0,0,3}; // >=    0-80%: 1-9
    const float centUp_Draw[ncent] = {9,7,6,4,2,7,4,9,7,5}; // <
    const float Nbin_Draw[nCent_Draw] = {938.80170, 579.89409, 288.35051, 91.37100, 21.37396, 386.08527, 56.99229, 301.05848,197, 168};
    const float NbinErr_Draw[nCent_Draw] = {26.28048, 28.86320, 30.39279, 21.05409, 8.93878, 29.86811, 14.85246, 23.32719, 1, 1};
    const float Npart_Draw[nCent_Draw] = {319.4, 227.6, 137.6, 60.5, 20.4, 167.82901, 40.83275, 126.49042, 1, 1};
    const float NpartErr_Draw[nCent_Draw] = {3.38968, 7.96962, 10.41932, 10.11154, 6.61248, 9.56374, 8.21022, 8.39142, 1, 1};

    
    const char namePrint[250] = "AuAu @ 200 GeV";
    const float ptCut = 0;
    const float fraction = 42.*1e3/0.0388 / 2. / 2.; // 42mb, B.R. = 0.0388, -1<y<1-->dy=2, D0+D0bar = 2
    
    float pt_mean[npt], pt_err[npt];
    for(int ipt=0; ipt<npt; ipt++) {
        pt_mean[ipt] = 0.5*(nptbin[ipt] + nptbin[ipt+1]);
        pt_err[ipt] = 0.5*(-nptbin[ipt] + nptbin[ipt+1]);
    }
    
    //read centrality
    // TFile* fin = new TFile("../D0_data_mix.root");
    TFile* fin = new TFile("../D0_Data_Mix_MergeUsed.root");
    TH1F* hcent = (TH1F*)fin->Get("hCentralityWeighted");
    hcent->SetDirectory(0);
    fin->Close();
    
    TFile* fPionEmd = new TFile("../pid_sys/piplus_sys.root");
    TF1* fpiTpc[nCent_Draw];
    TF1* fpiTpcRcp1[nCent_Draw];
    TF1* fpiTpcRcp2[nCent_Draw];
    for(int icent=0; icent<nCent_Draw; icent++) {
        fpiTpc[icent] = (TF1*)fPionEmd->Get(Form("fCombine3_%s",nameCent1_Draw[icent]));
        fpiTpc[icent] ->SetName(Form("fpiTpc_%s",nameCent1_Draw[icent]));
        fpiTpcRcp1[icent] = (TF1*)fPionEmd->Get(Form("fCombine3Rcp1_%s",nameCent1_Draw[icent]));
        fpiTpcRcp1[icent] ->SetName(Form("fpiTpcRcp1_%s",nameCent1_Draw[icent]));
        fpiTpcRcp2[icent] = (TF1*)fPionEmd->Get(Form("fCombine3Rcp2_%s",nameCent1_Draw[icent]));
        fpiTpcRcp2[icent] ->SetName(Form("fpiTpcRcp2_%s",nameCent1_Draw[icent]));
    }
    fPionEmd->Close();
    
    TFile* fKaonEmd = new TFile("../pid_sys/kaonminus_sys.root");
    TF1* fkTpc[nCent_Draw];
    TF1* fkTpcRcp1[nCent_Draw];
    TF1* fkTpcRcp2[nCent_Draw];
    for(int icent=0; icent<nCent_Draw; icent++) {
        fkTpc[icent] = (TF1*)fKaonEmd->Get(Form("fCombine3_%s",nameCent1_Draw[icent]));
        fkTpc[icent] ->SetName(Form("fkTpc_%s",nameCent1_Draw[icent]));
        fkTpcRcp1[icent] = (TF1*)fKaonEmd->Get(Form("fCombine3Rcp1_%s",nameCent1_Draw[icent]));
        fkTpcRcp1[icent] ->SetName(Form("fkTpcRcp1_%s",nameCent1_Draw[icent]));
        fkTpcRcp2[icent] = (TF1*)fKaonEmd->Get(Form("fCombine3Rcp2_%s",nameCent1_Draw[icent]));
        fkTpcRcp2[icent] ->SetName(Form("fkTpcRcp2_%s",nameCent1_Draw[icent]));
    }
    fKaonEmd->Close();
    
    //Convolution_D0_TPCsys.root
    TFile* fD0Emd = new TFile("../pid_sys/Convolution_D0_TPCsys.root");
    TProfile* fD0Tpc[ncent];
    TProfile* fD0TpcRcp1[ncent];
    TProfile* fD0TpcRcp2[ncent];
    for(int icent=0; icent<ncent; icent++) {
      fD0Tpc[icent] = (TProfile*)fD0Emd->Get(Form("h2Pt1Rebin_%s_pfx",nameCent1[icent]));
      fD0Tpc[icent] ->SetName(Form("fD0Tpc_%s",nameCent1[icent]));
      fD0TpcRcp1[icent] = (TProfile*)fD0Emd->Get(Form("h2Pt2Rebin_%s_pfx",nameCent1[icent]));
      fD0TpcRcp1[icent] ->SetName(Form("fD0TpcRcp1_%s",nameCent1[icent]));
      fD0TpcRcp2[icent] = (TProfile*)fD0Emd->Get(Form("h2Pt3Rebin_%s_pfx",nameCent1[icent]));
      fD0TpcRcp2[icent] ->SetName(Form("fD0TpcRcp2_%s",nameCent1[icent]));
      fD0Tpc[icent]->SetDirectory(0);
      fD0TpcRcp1[icent]->SetDirectory(0);
      fD0TpcRcp2[icent]->SetDirectory(0);
    }
    fD0Emd->Close();
    TFile* filepipid = new TFile("../pid_sys/pion_PidEff_Ks_170822_use_forSys.root");
    TF1* fpipid;
    fpipid = (TF1*)filepipid->Get("fpionNsigTof_effratio");
    fpipid ->SetName(Form("fpipid"));//9.87102992642156396e-01
    filepipid->Close();
    
    TFile* filekaonpid = new TFile("../pid_sys/kaon_PidEff_phi_170822_use_forSys.root");
    TF1* fkpid;
    fkpid = (TF1*)filekaonpid->Get("fkaonNsigTof_effratio");
    fkpid ->SetName(Form("fkpid")); // 9.87039958588842303e-01
    filekaonpid->Close();
    
    
    // calculate D0 cross section at pt > 0
    // Sys1: not correlated sys error between different pT range -- raw yield extraction
    // Sys2: correlated sys error between different pT range -- other sys error
    // Sys3: Nbinary error
    
    float XC[nCent_Draw], XC_Sys1[nCent_Draw], XC_Sys2[nCent_Draw], XC_Sys3[nCent_Draw], XC_Stat[nCent_Draw], XC_Sys[nCent_Draw]; //
    float XC_SysTest[nCent_Draw]; // test the last data piont's sys error
    // init
    for(int icent=0; icent<nCent_Draw; icent++) {
        XC[icent] = 0;
        XC_Sys[icent] = 0;
        XC_Sys1[icent] = 0;
        XC_Sys2[icent] = 0;
        XC_Sys3[icent] = 0;
        XC_Stat[icent] = 0;
        XC_SysTest[icent] = 0;
    }
    // calculte
    for(int icent=0; icent<nCent_Draw; icent++) {
        ifstream in;
        
        float yield[npt], yieldErr[npt], yieldSys[npt], yieldSys1[npt], yieldSys2[npt], yieldSys3[npt], tmp;
        float XC_SysTmp = 0;
        
        // mean value and stat error, double count sys error
        in.open(Form("../DoubleCount/data/yieldSys_%s.txt",nameCent1_Draw[icent]));
        if(in.eof()) { cout << "No double count sys error file!!!" << endl; exit(1);}
        XC_SysTmp = 0;
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            in >> yield[ipt] >> yieldErr[ipt] >> yieldSys[ipt] >> tmp;
            if(pt_mean[ipt]>ptCut) XC[icent] += yield[ipt];
            if(pt_mean[ipt]>ptCut) XC_Stat[icent] += pow(yieldErr[ipt],2);
            if(pt_mean[ipt]>ptCut) XC_SysTmp += yield[ipt]*yieldSys[ipt];
            if(ipt == npt-1) XC_SysTest[icent] += pow(yieldSys[ipt],2);
        }
        // XC_Sys2[icent] += pow(XC_SysTmp,2);
        // remobe double count for Rcp sys
        in.close();
        
        // sys 1 -- tpc track
        XC_SysTmp = 0;
        // double TmpCom =  sqrt(pow(fpiTpc[icent]->Eval(1),2) + pow(1.-fpipid->Eval(1),2)) + sqrt(pow(fkTpc[icent]->Eval(1),2) + pow(1.-fkpid->Eval(1),2));
        double TmpCom[npt];
        for(int ipt=0; ipt<npt; ipt++) {
          // TmpCom[ipt] =  sqrt(pow(fD0Tpc[icent]->GetBinContent(fD0Tpc[icent]->GetXaxis()->FindBin(pt_mean[ipt])),2) + pow(1.-fpipid->Eval(1),2) + pow(1.-fkpid->Eval(1),2));
          TmpCom[ipt] =  sqrt(pow(fD0TpcRcp1[icent]->GetBinContent(fD0TpcRcp1[icent]->GetXaxis()->FindBin(pt_mean[ipt])),2));// PID was cancal
          cout<< "TmpCom[ipt] = " << TmpCom[ipt] << endl;
        }
        
        for(int ipt=0; ipt<npt; ipt++) {
            yieldSys[ipt] = TmpCom[ipt];
            if(pt_mean[ipt]>ptCut) XC_SysTmp += yield[ipt]*yieldSys[ipt];
            if(ipt == npt-1) XC_SysTest[icent] += pow(yieldSys[ipt],2);
        }
        XC_Sys2[icent] += pow(XC_SysTmp,2);
        cout << "after TPC: "<<sqrt(XC_Sys2[icent])/XC[icent] << endl;
        
        //sys 2
        //-- 2.1 count with side band
        in.open(Form("../count/data/yieldSys_%s.txt",nameCent1_Draw[icent]));
        if(in.eof()) { cout << "No cout sys error file!!!" << endl; exit(1);}
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            in >> tmp >> yieldSys1[ipt];
        }
        in.close();
        
        //-- 2.2 count with fit change fit range
        in.open(Form("../fitRange/data/yieldSys_%s.txt",nameCent1_Draw[icent]));
        if(in.eof()) { cout << "No change fit range sys error file!!!" << endl; exit(1);}
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            in >> tmp >> yieldSys2[ipt];
        }
        in.close();
        
        //-- 2.3 count with likesign bkg subtraction
        in.open(Form("../likeSign/data/yieldSys_%s.txt",nameCent1_Draw[icent]));
        if(in.eof()) { cout << "No like-sign sys error file!!!" << endl; exit(1);}
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            in >> tmp >> yieldSys3[ipt];
        }
        in.close();
        
        for(int ipt=0; ipt<npt; ipt++) {
            yieldSys[ipt] = yieldSys1[ipt] > yieldSys2[ipt] ? yieldSys1[ipt] : yieldSys2[ipt];
            yieldSys[ipt] = yieldSys[ipt] > yieldSys3[ipt] ? yieldSys[ipt] : yieldSys3[ipt];
            if(pt_mean[ipt]>ptCut) XC_Sys1[icent] += pow(yieldSys[ipt]*yield[ipt],2);
            if(ipt == npt-1) XC_SysTest[icent] += pow(yieldSys[ipt],2);
        }
        cout << "Yield : "<<sqrt(XC_Sys1[icent])/XC[icent] << endl;
        
        //sys 3, Daughter pt Cut scan // choose the maximum difference
        //sys 3.1 -- daughter pt cut1
        in.open(Form("../ptCut1/data/yieldSys_%s.txt",nameCent1_Draw[icent]));
        if(in.eof()) { cout << "No pt = 0.3 sys error file!!!" << endl; exit(1);}
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            in >> tmp >> yieldSys1[ipt];
        }
        in.close();
        
        //sys 3.2 -- daughter pt cut2
        in.open(Form("../ptCut2/data/yieldSys_%s.txt",nameCent1_Draw[icent]));
        if(in.eof()) { cout << "No pt = 0.5 sys error file!!!" << endl; exit(1);}
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            in >> tmp >> yieldSys2[ipt];
        }
        in.close();
        
        XC_SysTmp = 0;
        for(int ipt=0; ipt<npt; ipt++) {
            yieldSys[ipt] = yieldSys1[ipt] > yieldSys2[ipt] ? yieldSys1[ipt] : yieldSys2[ipt];
            if(pt_mean[ipt]>ptCut) XC_SysTmp += yield[ipt]*yieldSys[ipt];
            if(ipt == npt-1) XC_SysTest[icent] += pow(yieldSys[ipt],2);
        }
        // XC_Sys2[icent] += pow(XC_SysTmp,2);// correlated removed to CrossX_Rcp_sys
        
        //sys 4 topological cut
        //4.1 -- tight topo cuts
        in.open(Form("../topoCut1/data/yieldSys_%s.txt",nameCent1_Draw[icent]));
        if(in.eof()) { cout << "No tight topo sys error file!!!" << endl; exit(1);}
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            in >> tmp >> yieldSys1[ipt];
        }
        in.close();
        
        //4.2 -- loose topo cuts
        in.open(Form("../topoCut2/data/yieldSys_%s.txt",nameCent1_Draw[icent]));
        if(in.eof()) { cout << "No loose topo sys error file!!!" << endl; exit(1);}
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            in >> tmp >> yieldSys2[ipt];
        }
        in.close();
        
        XC_SysTmp = 0;
        for(int ipt=0; ipt<npt; ipt++) {
            // tmp[ipt] = tmp1[ipt] > tmp2[ipt] ? tmp1[ipt] : tmp2[ipt];
            yieldSys[ipt] = (yieldSys1[ipt] + yieldSys2[ipt])/2.0;;
            if(pt_mean[ipt]>ptCut) XC_SysTmp += yield[ipt]*yieldSys[ipt];
            if(ipt == npt-1) XC_SysTest[icent] += pow(yieldSys[ipt],2);
        }
        // XC_Sys2[icent] += pow(XC_SysTmp,2);// correlated removed to CrossX_Rcp_sys
        
        //5 do the vertex reso. correction
        // vtx sys. error
        in.open(Form("../vtxCorr/data/vtxSys_%s.txt",nameCent1_Draw[icent]));
        if(in.eof()) { cout << "No vtx sys error file!!!" << endl; exit(1);}
        XC_SysTmp = 0;
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            in >> yieldSys[ipt];
            if(pt_mean[ipt]>ptCut) XC_SysTmp += yield[ipt]*yieldSys[ipt];
            if(ipt == npt-1) XC_SysTest[icent] += pow(yieldSys[ipt],2);
        }
        XC_Sys2[icent] += pow(XC_SysTmp,2);
        in.close();
        
        // vtx stat. error
        in.open(Form("../vtxCorr/data/vtxStat_%s.txt",nameCent1_Draw[icent]));
        if(in.eof()) { cout << "No vtx sys error file!!!" << endl; exit(1);}
        XC_SysTmp = 0;
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            in >> yieldSys[ipt];
            if(pt_mean[ipt]>ptCut) XC_SysTmp += yield[ipt]*yieldSys[ipt];
            if(ipt == npt-1) XC_SysTest[icent] += pow(yieldSys[ipt],2);
        }
        XC_Sys2[icent] += pow(XC_SysTmp,2);
        in.close();

        cout << "after vtx : "<<sqrt(XC_Sys2[icent])/XC[icent] << endl;
        
        //6 Nbin error
        XC_Sys3[icent] = XC[icent]*NbinErr_Draw[icent]/Nbin_Draw[icent];
        
        
        // convert to cross section: mub
        int bin_lw = hcent->FindBin(centLw_Draw[icent]+1e-6);
        int bin_up = hcent->FindBin(centUp_Draw[icent]-1e-6);
        //cout << "cent bin: " << bin_lw << "-" << bin_up << endl;
        float nevents = hcent->Integral(bin_lw,bin_up);
        XC[icent] = XC[icent] * fraction / nevents / Nbin_Draw[icent];
        XC_Stat[icent] = sqrt(XC_Stat[icent]) * fraction / nevents / Nbin_Draw[icent];
        XC_Sys1[icent] = sqrt(XC_Sys1[icent]) * fraction / nevents / Nbin_Draw[icent];
        XC_Sys2[icent] = sqrt(XC_Sys2[icent]) * fraction / nevents / Nbin_Draw[icent];
        XC_Sys3[icent] = XC_Sys3[icent] * fraction / nevents / Nbin_Draw[icent];
        
        XC_Sys[icent] = sqrt(pow(XC_Sys1[icent],2) + pow(XC_Sys2[icent],2));
        XC_SysTest[icent] = sqrt(XC_SysTest[icent]);
        
        cout << nameCent_Draw[icent] << endl;
        cout << "Mean = " << XC[icent] << ", stat.err(%) = " << XC_Stat[icent]/XC[icent] << ", sysNotCorr.err(%) = " << XC_Sys1[icent]/XC[icent] << ", sysCorr.err(%) = " << XC_Sys2[icent]/XC[icent] << ", sys.err(%) = " << XC_Sys[icent]/XC[icent] << ", sysNbin.err(%) = " << XC_Sys3[icent]/XC[icent] << endl;
        cout << "last pT bin's sys error(%): " << XC_SysTest[icent] << endl;
    }

    ofstream outData;
    outData.open(Form("D0_Xsection_forRcp.txt"));
    outData << "centrality, " << " Nbin, " << " Xsection, " << " stat. " << "sys."<< " sys1.(NotCorr), " << " sys2(Corr). " << "sys3(Nbin)."<< endl;
    for(int icent=0; icent<nCent_Draw; icent++) {
      outData << nameCent_Draw[icent] << " " << Nbin_Draw[icent] << " "<< XC[icent] << " " << XC_Stat[icent] <<" " << XC_Sys[icent] <<" "  << XC_Sys1[icent] << " " <<XC_Sys2[icent] <<" " << XC_Sys3[icent] << endl;
    }
    outData.close();


    // TGraphErrors* gD0CorssS_err = new TGraphErrors(nCent_Draw, Nbin_Draw, XC, 0, XC_Stat);
    // TGraphErrors* gD0CorssS_sys = new TGraphErrors(nCent_Draw, Nbin_Draw, XC, 0, XC_Sys);
    // TGraphErrors* gD0CorssS_sys1 = new TGraphErrors(nCent_Draw, Nbin_Draw, XC, 0, XC_Sys3);
    //
    // // write
    // sprintf(dir,"data");
    // sprintf(CMD,"[ -d %s ] || mkdir -p %s",dir,dir);
    // gSystem->Exec(CMD);
    //
    // sprintf(name,"%s/D0_CrossSection_ptGT%d.root",dir,ptCut);
    // TFile* fout = new TFile(name,"RECREATE");
    // sprintf(name,"D0CrossSection_err_ptGT%d",ptCut);
    // gD0CorssS_err->Write(name);
    // sprintf(name,"D0CrossSection_sys_ptGT%d",ptCut);
    // gD0CorssS_sys->Write(name);
    // sprintf(name,"D0CrossSection_NbinError_ptGT%d",ptCut);
    // gD0CorssS_sys1->Write(name);
    // fout->Close();
    //
    // // plot
    // sprintf(dir,"pic");
    // sprintf(CMD,"[ -d %s ] || mkdir -p %s",dir,dir);
    // gSystem->Exec(CMD);
    //
    // TCanvas* c1 = new TCanvas("c1", "A Canvas",10,10,1000,1000);
    // setPad(c1);
    //
    // // find ymin/ymax
    // float ymin = 9e6;
    // float ymax = -9e6;
    // for(int icent=0; icent<nCent_Draw; icent++) {
    //     if(ymax<XC[icent]) ymax = XC[icent];
    //     if(ymin>XC[icent]) ymin = XC[icent];
    // }
    //
    // gPad->SetLogx();
    // h0 = new TH1F("","",1,5,Nbin_Draw[nCent_Draw-1]*5.);
    // h0->GetYaxis()->SetRangeUser(0,ymax*1.5);
    // h0->SetLineColor(kWhite);
    // setHisto(h0,"","N_{bin}", "d#sigma/dy|_{y=0} (#mub)");
    // h0->Draw();
    //
    // sprintf(name,"%s",namePrint);
    // drawLatex(0.62,0.88,name,132,0.045,1);
    // sprintf(name,"p_{T} > %.0f GeV/c",ptCut);
    // drawLatex(0.62,0.82,name,132,0.045,1);
    //
    // const float lineWidth = 3;
    // const float markersize = 2.0;
    //
    // // draw Nbin error
    // const float boxW = 0.15;
    // for(int i=0; i<gD0CorssS_sys1->GetN(); i++) {
    //     TBox* tb = new TBox(gD0CorssS_sys1->GetX()[i]*(1-boxW),gD0CorssS_sys1->GetY()[i]-gD0CorssS_sys1->GetEY()[i],gD0CorssS_sys1->GetX()[i]*(1+boxW)*1.0,gD0CorssS_sys1->GetY()[i]+gD0CorssS_sys1->GetEY()[i]);
    //     tb->SetFillColor(kGreen);
    //     tb->SetLineColor(kGreen);
    //     tb->Draw("same");
    //     if(!legend) {
    //         legend = new TLegend(0.18,0.14,0.4,0.18);
    //         legend->SetFillStyle(0);
    //         legend->SetFillColor(10);
    //         legend->SetBorderSize(0);
    //         legend->SetTextSize(0.03);
    //         legend->SetTextFont(132);
    //         legend->AddEntry(tb,"N_{bin} error","f");
    //         legend->Draw();
    //     }
    // }
    //
    //
    // //draw systematic error
    // const float sysw = 0.1;
    // for(int i=0; i<gD0CorssS_sys->GetN(); i++) {
    //     const float sysl = gD0CorssS_sys->GetY()[i] * 0.15;
    //     TLine *llw = new TLine(gD0CorssS_sys->GetX()[i]*(1-sysw),gD0CorssS_sys->GetY()[i]-gD0CorssS_sys->GetEY()[i],gD0CorssS_sys->GetX()[i]*(1+sysw),gD0CorssS_sys->GetY()[i]-gD0CorssS_sys->GetEY()[i]);
    //     llw->SetLineWidth(lineWidth);
    //     llw->SetLineColor(kBlack);
    //     llw->Draw("same");
    //     TLine *lhi = new TLine(gD0CorssS_sys->GetX()[i]*(1-sysw),gD0CorssS_sys->GetY()[i]+gD0CorssS_sys->GetEY()[i],gD0CorssS_sys->GetX()[i]*(1+sysw),gD0CorssS_sys->GetY()[i]+gD0CorssS_sys->GetEY()[i]);
    //     lhi->SetLineWidth(lineWidth);
    //     lhi->SetLineColor(kBlack);
    //     lhi->Draw("same");
    //     TLine *lv1 = new TLine(gD0CorssS_sys->GetX()[i]*(1-sysw),gD0CorssS_sys->GetY()[i]-gD0CorssS_sys->GetEY()[i],gD0CorssS_sys->GetX()[i]*(1-sysw),gD0CorssS_sys->GetY()[i]+gD0CorssS_sys->GetEY()[i]);
    //     lv1->SetLineWidth(lineWidth);
    //     lv1->SetLineColor(kBlack);
    //     lv1->Draw("same");
    //     TLine *lv2 = new TLine(gD0CorssS_sys->GetX()[i]*(1+sysw),gD0CorssS_sys->GetY()[i]-gD0CorssS_sys->GetEY()[i],gD0CorssS_sys->GetX()[i]*(1+sysw),gD0CorssS_sys->GetY()[i]+gD0CorssS_sys->GetEY()[i]);
    //     lv2->SetLineWidth(lineWidth);
    //     lv2->SetLineColor(kBlack);
    //     lv2->Draw("same");
    // }
    //
    // gD0CorssS_err->SetMarkerColor(kBlack);
    // gD0CorssS_err->SetMarkerStyle(kFullCircle);
    // gD0CorssS_err->SetMarkerSize(markersize);
    // gD0CorssS_err->SetLineColor(kBlack);
    // gD0CorssS_err->SetLineWidth(lineWidth);
    // gD0CorssS_err->Draw("psame");
    //
    // sprintf(name,"%s/D0_CrossS_ptGT%d.gif",dir,ptCut);
    // c1->SaveAs(name);
}
