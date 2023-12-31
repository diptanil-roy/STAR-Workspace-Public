#include "../myFunction.h"
#include "../myConst.h"
#include "../anaCuts.h"
void plot_err() {
    globalSetting();
    char dir[250];
    char name[250];
    char title[250];
    //TString CMD = 0;
    char CMD[250];
    TLegend* legend;
    TH1F* h0;
    
    sprintf(dir,"pic");
    sprintf(CMD,"[ -d %s ] || mkdir -p %s",dir,dir);
    gSystem->Exec(CMD);
    
    float ptmean[npt], pterr[npt];
    for(int ipt=0; ipt<npt; ipt++) {
        ptmean[ipt] = 0.5*(nptbin[ipt]+nptbin[ipt+1]);
        pterr[ipt] = 0.5*(-nptbin[ipt]+nptbin[ipt+1]);
    }
    
    // const int nerr = 6;
    // const char nameErr[nerr][250] = {"TPC track", "raw yield extract", "single track p_{T}", "50% topo. eff.", "150% topo. eff.", "double count"};
    //
    // const int nerr = 7;
    // const char nameErr[nerr][250] = {"TPC track", "raw yield extract", "single track p_{T}1","single track p_{T}2", "50% topo. eff.", "150% topo. eff.", "double count"};
    
    // const int nerr = 6;
    // const char nameErr[nerr][250] = {"TPC track", "raw yield extract", "single track p_{T}", "50%/150% topo. eff.", "double count", "vertex correction"};

    // const int nerr = 7;
    // const char nameErr[nerr][250] = {"TPC track", "PID", "raw yield extract", "single track p_{T}", "50%/150% topo. eff.", "double count", "vertex correction"};

    // const int nerr = 8;
    // const char nameErr[nerr][250] = {"TPC track", "PID", "raw yield extract", "single track p_{T}", "50%/150% topo. eff.", "double count", "vertex correction", "secondaryTrack"};
    //
    const int nerr = 10;
    const char nameErr[nerr][250] = {"TPC track", "PID", "raw yield extract", "single track p_{T}", "50%/150% topo. eff.", "double count", "vertex correction", "secondaryTrack", "tpc luminosity","general sys."};
    
    float pt_mean[npt], pt_err[npt];
    for(int i=0; i<npt; i++) {
        pt_mean[i] = 0.5*(nptbin[i] + nptbin[i+1]);
        pt_err[i] = 0.5*(-nptbin[i] + nptbin[i+1]);
    }
    
    TFile* fPionEmd = new TFile("../pid_sys/piplus_sys.root");
    TF1* fpiTpc[ncent];
    TF1* fpiTpcRcp1[ncent];
    TF1* fpiTpcRcp2[ncent];
    for(int icent=0; icent<ncent; icent++) {
      fpiTpc[icent] = (TF1*)fPionEmd->Get(Form("fCombine3_%s",nameCent1[icent]));
      fpiTpc[icent] ->SetName(Form("fpiTpc_%s",nameCent1[icent]));
      fpiTpcRcp1[icent] = (TF1*)fPionEmd->Get(Form("fCombine3Rcp1_%s",nameCent1[icent]));
      fpiTpcRcp1[icent] ->SetName(Form("fpiTpcRcp1_%s",nameCent1[icent]));
      fpiTpcRcp2[icent] = (TF1*)fPionEmd->Get(Form("fCombine3Rcp2_%s",nameCent1[icent]));
      fpiTpcRcp2[icent] ->SetName(Form("fpiTpcRcp2_%s",nameCent1[icent]));
    }
    fPionEmd->Close();

    TFile* fKaonEmd = new TFile("../pid_sys/kaonminus_sys.root");
    TF1* fkTpc[ncent];
    TF1* fkTpcRcp1[ncent];
    TF1* fkTpcRcp2[ncent];
    for(int icent=0; icent<ncent; icent++) {
      fkTpc[icent] = (TF1*)fKaonEmd->Get(Form("fCombine3_%s",nameCent1[icent]));
      fkTpc[icent] ->SetName(Form("fkTpc_%s",nameCent1[icent]));
      fkTpcRcp1[icent] = (TF1*)fKaonEmd->Get(Form("fCombine3Rcp1_%s",nameCent1[icent]));
      fkTpcRcp1[icent] ->SetName(Form("fkTpcRcp1_%s",nameCent1[icent]));
      fkTpcRcp2[icent] = (TF1*)fKaonEmd->Get(Form("fCombine3Rcp2_%s",nameCent1[icent]));
      fkTpcRcp2[icent] ->SetName(Form("fkTpcRcp2_%s",nameCent1[icent]));
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

    // read sys err
    ifstream in;
    float systmp, ytmp;
    float sys[ncent][nerr][npt];
    float y[npt], yerr[npt], ysys[npt];
    float tmp[npt];
    float tmp1[npt];
    float tmp2[npt];
    float tmp3[npt];
    for(int icent=0; icent<ncent; icent++) {
        int ierr = 0;
        
        // double TmpCom =  sqrt(pow(fpiTpc[icent]->Eval(1),2) + pow(1.-fpipid->Eval(1),2)) + sqrt(pow(fkTpc[icent]->Eval(1),2) + pow(1.-fkpid->Eval(1),2));
        // double TmpCom =  sqrt(pow(fpiTpc[icent]->Eval(1),2)) + sqrt(pow(fkTpc[icent]->Eval(1),2));//seperate PID
        double TmpCom[npt];
        for(int ipt=0; ipt<npt; ipt++) {
          TmpCom[ipt] =  fD0Tpc[icent]->GetBinContent(fD0Tpc[icent]->GetXaxis()->FindBin(ptmean[ipt]));
          cout<< "TmpCom[ipt] = " << TmpCom[ipt] << endl;
        }

        //sys 1 -- tpc track
        for(int ipt=0; ipt<npt; ipt++) {
            sys[icent][ierr][ipt] = TmpCom[ipt];
            // sys[icent][ierr][ipt] = 0.1;
            //ysys[ipt] = sqrt(pow(0.04,2)+pow(ysys[ipt],2));
        }
        ierr++;
        
        // double TmpCom =  sqrt(pow(fpiTpc[icent]->Eval(1),2) + pow(1.-fpipid->Eval(1),2)) + sqrt(pow(fkTpc[icent]->Eval(1),2) + pow(1.-fkpid->Eval(1),2));
        double TmpCom2 =  sqrt(pow(1.-fpipid->Eval(1),2)) + sqrt(pow(1.-fkpid->Eval(1),2));
        //sys 1.2 -- pid seperated
        for(int ipt=0; ipt<npt; ipt++) {
            sys[icent][ierr][ipt] = TmpCom2;
        }
        ierr++;
        
        //sys 2 
        //-- 2.1 count with side band
        in.open(Form("../count/data/yieldSys_%s.txt",nameCent1[icent]));
        if(in.eof()) { cout << "No cout sys error file!!!" << endl; exit(1);}
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            in >> y[ipt] >> tmp1[ipt];
        }
        in.close();
        
        //-- 2.2 count with fit change fit range
        in.open(Form("../fitRange/data/yieldSys_%s.txt",nameCent1[icent]));
        if(in.eof()) { cout << "No change fit range sys error file!!!" << endl; exit(1);}
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            in >> y[ipt] >> tmp2[ipt];
        }
        in.close();

        //-- 2.3 count with likesign bkg subtraction
        in.open(Form("../likeSign/data/yieldSys_%s.txt",nameCent1[icent]));
        if(in.eof()) { cout << "No like-sign sys error file!!!" << endl; exit(1);}
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            in >> y[ipt] >> tmp3[ipt];
        }
        in.close();

        for(int ipt=0; ipt<npt; ipt++) {
          tmp[ipt] = tmp1[ipt] > tmp2[ipt] ? tmp1[ipt] : tmp2[ipt];
          tmp[ipt] = tmp[ipt] > tmp3[ipt] ? tmp[ipt] : tmp3[ipt];
          // ysys[ipt] = sqrt(pow(tmp[ipt],2)+pow(ysys[ipt],2));
        }

        for(int ipt=0; ipt<npt; ipt++) {
          sys[icent][ierr][ipt] = tmp[ipt];
        }
        ierr++;

        //sys 3, Daughter pt Cut scan // choose the maximum difference
        //sys 3.1 -- daughter pt cut1
        in.open(Form("../ptCut1/data/yieldSys_%s.txt",nameCent1[icent]));
        if(in.eof()) { cout << "No pt = 0.3 sys error file!!!" << endl; exit(1);}
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            in >> y[ipt] >> tmp1[ipt];
        }
        in.close();
        
        //sys 3.2 -- daughter pt cut2
        in.open(Form("../ptCut2/data/yieldSys_%s.txt",nameCent1[icent]));
        if(in.eof()) { cout << "No pt = 0.5 sys error file!!!" << endl; exit(1);}
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            in >> y[ipt] >> tmp2[ipt];
        }
        in.close();

        for(int ipt=0; ipt<npt; ipt++) {
          tmp[ipt] = tmp1[ipt] > tmp2[ipt] ? tmp1[ipt] : tmp2[ipt];
          // ysys[ipt] = sqrt(pow(tmp[ipt],2)+pow(ysys[ipt],2));
        }

        for(int ipt=0; ipt<npt; ipt++) {
          sys[icent][ierr][ipt] = tmp[ipt];
        }
        ierr++;

        
        //sys 4 topological cut
        //4.1 -- tight topo cuts
        in.open(Form("../topoCut1/data/yieldSys_%s.txt",nameCent1[icent]));
        if(in.eof()) { cout << "No tight topo sys error file!!!" << endl; exit(1);}
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            in >> y[ipt] >> tmp1[ipt];
            //if(icent == 4 && ipt >= npt-3) tmp1[ipt] =0; // for the tight cut, the signal is not good, no need to include this sys for some certain pt bin
        }
        in.close();
        
        //4.2 -- loose topo cuts
        in.open(Form("../topoCut2/data/yieldSys_%s.txt",nameCent1[icent]));
        if(in.eof()) { cout << "No loose topo sys error file!!!" << endl; exit(1);}
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            in >> y[ipt] >> tmp2[ipt];
        }
        in.close();
        
        for(int ipt=0; ipt<npt; ipt++) {
           tmp[ipt] = tmp1[ipt] > tmp2[ipt] ? tmp1[ipt] : tmp2[ipt];
          // ysys[ipt] = sqrt(pow(tmp[ipt],2)+pow(ysys[ipt],2));
        }
        
        for(int ipt=0; ipt<npt; ipt++) {
           // sys[icent][ierr][ipt] = tmp[ipt];
           sys[icent][ierr][ipt] = (tmp1[ipt]+tmp2[ipt])/2.0;
        }
        ierr++;

        //sys 5 -- double count
        in.open(Form("../DoubleCount/data/yieldSys_%s.txt",nameCent1[icent]));
        if(in.eof()) { cout << "No double count sys error file!!!" << endl; exit(1);}
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            float errtmp;
            float tmptmp;
            in >> ytmp >> errtmp >> systmp >> tmptmp;
            sys[icent][ierr][ipt] = systmp;
        }
        in.close();
        ierr++;
        
        //sys 6 vtx resolution
        //4.1 -- vtx sys.
        in.open(Form("../vtxCorr/data/vtxSys_%s.txt",nameCent1[icent]));
        if(in.eof()) { cout << "No vtx sys error file!!!" << endl; exit(1);}
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            in >> tmp1[ipt];
        }
        in.close();
        
        //4.2 -- vtx stat.
        in.open(Form("../vtxCorr/data/vtxStat_%s.txt",nameCent1[icent]));
        if(in.eof()) { cout << "No vtx sys error file!!!" << endl; exit(1);}
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            in >> tmp2[ipt];
        }
        in.close();
        
        for(int ipt=0; ipt<npt; ipt++) {
            sys[icent][ierr][ipt] = sqrt(pow(tmp1[ipt],2)+pow(tmp2[ipt],2) - 0.05*0.05);// 5% general sys was subtracted here
        }
        ierr++;
        
        //sys 7 second track
        // for(int ipt=0; ipt<npt; ipt++) {
            // sys[icent][ierr][ipt] = 0.02;
        // }
        // ierr++;

        // secondary sys from the pions 
        in.open(Form("../secondarySys/data/yieldSys_%s.txt",nameCent1[icent]));
        if(in.eof()) { cout << "No vtx sys error file!!!" << endl; exit(1);}
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            in >> tmp[ipt] >> tmp[ipt];
            sys[icent][ierr][ipt] = sqrt(pow(tmp[ipt],2)+pow(0.02,2));
        }
        in.close();
        ierr++;
        
        //sys  tpc luminosity dep
        for(int ipt=0; ipt<npt; ipt++) {
            sys[icent][ierr][ipt] = 0.04;
        }
        ierr++;

        //sys 8 second track
        for(int ipt=0; ipt<npt; ipt++) {
            sys[icent][ierr][ipt] = 0.05;
        }
        ierr++;
    }
    
    //define sys err graph
    TGraph* gSys[ncent][nerr];
    for(int icent=0; icent<ncent; icent++) {
        for(int ierr=0; ierr<nerr; ierr++) {
            gSys[icent][ierr] = new TGraph(npt,pt_mean,sys[icent][ierr]);
        }
    }

    for(int icent=0; icent<ncent; icent++) {
      if(!(icent == 0 || icent==4) ) continue; // print 0-10%

        for(int ierr=0; ierr<nerr; ierr++) {
          for(int ipt=0; ipt<npt; ipt++) {
            cout << nameCent[icent] << " " << nameErr[ierr] << " "<<  sys[icent][ierr][ipt] << endl;
          }
        }

          for(int ipt=0; ipt<npt; ipt++) {
            cout << " combine HFT " << nameCent[icent] << " "<<  sqrt(pow(sys[icent][3][ipt],2)+pow(sys[icent][4][ipt],2)) << endl;
          }

    }
    
    //plot
    TCanvas* c1 = new TCanvas("c1", "A Canvas",10,10,600,800);
    setPad(c1);
    gPad->SetLogy();
    for(int icent=0; icent<ncent; icent++) {
        float ymin = 1.01e-3;
        float ymax = 9.5e0;
        h0= new TH1F("","",1,0,nptbin[npt]);
        h0->Draw();
        h0->SetMinimum(ymin),
        h0->SetMaximum(ymax);
        setHisto(h0,"","p_{T} (GeV/c)", "Sys. error (#times 100%)");
        const float legy_lw = 0.93 - 0.04*nerr;
        legend = new TLegend(0.6,legy_lw,0.9,0.93);
        legend->SetFillStyle(0);
        legend->SetFillColor(10);
        legend->SetBorderSize(0);
        legend->SetTextSize(0.04);
        legend->SetTextFont(132);
        for(int ierr=0; ierr<nerr; ierr++) {
            legend->AddEntry(gSys[icent][ierr],nameErr[ierr],"p");
            gSys[icent][ierr]->Draw("psame");
            gSys[icent][ierr]->SetMarkerStyle(MARKERSTYLE[ierr]);
            gSys[icent][ierr]->SetMarkerColor(COLOR[ierr]);
            gSys[icent][ierr]->SetMarkerSize(1.5);
        }
        legend->Draw();
        drawLatex(0.18,0.9,Form("AuAu @200GeV  %s",nameCent[icent]),132,0.04,1);
        gPad->SetLogy();
        sprintf(name,"%s/sysErr_%s_2.pdf",dir,nameCent1[icent]);
        c1->SaveAs(name);
    }
}
