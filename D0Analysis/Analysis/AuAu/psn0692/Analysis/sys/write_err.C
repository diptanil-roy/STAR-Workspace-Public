#include "../anaCuts.h"
#include "../myConst.h"
void write_err() {
    char CMD[250];
    char dir[250];
    sprintf(dir,"data");
    sprintf(CMD,"[ -d %s ] || mkdir -p %s",dir,dir);
    gSystem->Exec(CMD);
    
    float ptmean[npt], pterr[npt];
    for(int ipt=0; ipt<npt; ipt++) {
        ptmean[ipt] = 0.5*(nptbin[ipt]+nptbin[ipt+1]);
        pterr[ipt] = 0.5*(-nptbin[ipt]+nptbin[ipt+1]);
    }
    
    //read centrality
    // TFile* fin = new TFile("../D0_data_mix.root");
    TFile* fin = new TFile("../D0_Data_Mix_MergeUsed.root");
    TH1F* hcent = (TH1F*)fin->Get("hCentralityWeighted");
    hcent->SetDirectory(0);
    fin->Close();

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

    //caculate err
    float y[npt], yerr[npt], ysys[npt];
    float tmp[npt];
    float tmp1[npt];
    float tmp2[npt];
    float tmp3[npt];
    ifstream in;
    ofstream out;
    TFile* fout = new TFile("D0_Spectra_Run14_HFT_beforePtShift.root","RECREATE");
    //TFile* fin = new TFile(
    for(int icent=0; icent<ncent; icent++) {
        //init
        for(int ipt=0; ipt<npt; ipt++) {
            ysys[ipt] = 0;
            yerr[ipt] = 0;
        }
        
        //statistics error
        // in.open(Form("../default/data/yield_%s.txt",nameCent1[icent]));
        in.open(Form("../default/data/re_yield_%s.txt",nameCent1[icent]));
        if(in.eof()) { cout << "No default yield error file!!!" << endl; exit(1);}
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) cout << {
                break;
            }
            in >> y[ipt] >> yerr[ipt];
            yerr[ipt] /= y[ipt];
        }
        in.close();

        // double TmpCom =  sqrt(pow(fpiTpc[icent]->Eval(1),2) + pow(1.-fpipid->Eval(1),2)) + sqrt(pow(fkTpc[icent]->Eval(1),2) + pow(1.-fkpid->Eval(1),2));
        double TmpCom[npt];
        for(int ipt=0; ipt<npt; ipt++) {
          TmpCom[ipt] =  sqrt(pow(fD0Tpc[icent]->GetBinContent(fD0Tpc[icent]->GetXaxis()->FindBin(ptmean[ipt])),2) + pow(1.-fpipid->Eval(1),2) + pow(1.-fkpid->Eval(1),2));
          cout<< "TmpCom[ipt] = " << TmpCom[ipt] << endl;
        }
        
        //sys 1 -- tpc track
        for(int ipt=0; ipt<npt; ipt++) {
            ysys[ipt] = sqrt(pow(TmpCom[ipt],2)+pow(ysys[ipt],2)+pow(0.014,2));
            // ysys[ipt] = sqrt(pow(0.1,2)+pow(ysys[ipt],2));
            //ysys[ipt] = sqrt(pow(0.04,2)+pow(ysys[ipt],2));
        }
        
        //sys 2 
        //-- 2.1 count with side band
        in.open(Form("../count/data/yieldSys_%s.txt",nameCent1[icent]));
        if(in.eof()) { cout << "No cout sys error file!!!" << endl; exit(1);}
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            in >> y[ipt] >> tmp1[ipt];
            // ysys[ipt] = sqrt(pow(tmp,2)+pow(ysys[ipt],2));
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
          ysys[ipt] = sqrt(pow(tmp[ipt],2)+pow(ysys[ipt],2));
        }

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
          ysys[ipt] = sqrt(pow(tmp[ipt],2)+pow(ysys[ipt],2));
        }
        
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
          // tmp[ipt] = tmp1[ipt] > tmp2[ipt] ? tmp1[ipt] : tmp2[ipt];
            tmp[ipt] = (tmp1[ipt] + tmp2[ipt])/2.0;;
          ysys[ipt] = sqrt(pow(tmp[ipt],2)+pow(ysys[ipt],2));
        }
        
        //sys 5 -- double count
        in.open(Form("../DoubleCount/data/yieldSys_%s.txt",nameCent1[icent]));
        if(in.eof()) { cout << "No double count sys error file!!!" << endl; exit(1);}
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            in >> y[ipt] >> yerr[ipt] >> tmp[ipt] >> tmp1[ipt];
            ysys[ipt] = sqrt(pow(tmp[ipt],2)+pow(ysys[ipt],2));
            yerr[ipt] /= y[ipt];
        }
        in.close();
        
        // do the vertex reso. correction
        // vtx sys. error
        in.open(Form("../vtxCorr/data/vtxSys_%s.txt",nameCent1[icent]));
        if(in.eof()) { cout << "No vtx sys error file!!!" << endl; exit(1);}
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            in >> tmp[ipt];
            ysys[ipt] = sqrt(pow(tmp[ipt],2)+pow(ysys[ipt],2));
        }
        in.close();
        // vtx stat. error
        in.open(Form("../vtxCorr/data/vtxStat_%s.txt",nameCent1[icent]));
        if(in.eof()) { cout << "No vtx sys error file!!!" << endl; exit(1);}
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            in >> tmp[ipt];
            ysys[ipt] = sqrt(pow(tmp[ipt],2)+pow(ysys[ipt],2));
            // yerr[ipt] = sqrt(pow(tmp[ipt],2)+pow(yerr[ipt],2));
        }
        in.close();
        
        // secondary sys 2%
        for(int ipt=0; ipt<npt; ipt++) {
            ysys[ipt] = sqrt(pow(0.02,2)+pow(ysys[ipt],2));
        }

        // luminosity sys 4%
        for(int ipt=0; ipt<npt; ipt++) {
            ysys[ipt] = sqrt(pow(0.04,2)+pow(ysys[ipt],2));
        }

        // secondary sys from the pions 
        in.open(Form("../secondarySys/data/yieldSys_%s.txt",nameCent1[icent]));
        if(in.eof()) { cout << "No vtx sys error file!!!" << endl; exit(1);}
        for(int ipt=0; ipt<npt; ipt++) {
            if(in.eof()) break;
            in >> tmp[ipt] >> tmp[ipt];
            ysys[ipt] = sqrt(pow(tmp[ipt],2)+pow(ysys[ipt],2));
        }
        in.close();
        
        //transfer to the dN/Nev2piPtdptdy
        int bin_lw = hcent->FindBin(centLw[icent]+1e-6);
        int bin_up = hcent->FindBin(centUp[icent]-1e-6);
        cout << "cent bin: " << bin_lw << "-" << bin_up << endl;
        float nevents = hcent->Integral(bin_lw,bin_up);
        cout << "number of events: " << nevents << endl;
        float factor = 1./ 2.; //delta y
        float factor = factor*1./2.; // (D0+D0bar)/2
        cout << factor << endl;
        factor /= nevents;
        factor = factor/0.0388; //branch ratio
        //factor = factor*1000.*42.;
        //factor /= NbinMean[icent];
        for(int ipt=0; ipt<npt; ipt++) {
            float pt = 0.5*(nptbin[ipt]+nptbin[ipt+1]);
            float ptWidth = nptbin[ipt+1] - nptbin[ipt];
            //cout << pt << "\t" << ptWidth << endl;
            float scale = factor/(2.*TMath::Pi()*pt*ptWidth);
            //float scale = factor;
            
            y[ipt] *= scale;
            yerr[ipt] *= y[ipt];
            ysys[ipt] *= y[ipt];
        }
        
        //out put
        out.open(Form("data/sys_%s.txt",nameCent1[icent]));
        out << "yield" << "\t\t" << "sta.err" << "\t\t" << "sys.err" << "\t\t" << "sta.relative" << "\t" << "sys.relative" << endl;
        for(int ipt=0; ipt<npt; ipt++) {
            out << y[ipt] << "\t" << yerr[ipt] << "\t" << ysys[ipt] << "\t";
            out << yerr[ipt]/y[ipt] << "\t" << ysys[ipt]/y[ipt] << endl;
        }
        out.close();
        
        fout->cd();
        TGraphErrors* gD0err = new TGraphErrors(npt, ptmean, y, 0, yerr);
        TGraphErrors* gD0sys = new TGraphErrors(npt, ptmean, y, 0, ysys);
        gD0err->SetMarkerStyle(MARKERSTYLE[0]);
        gD0err->SetMarkerSize(2.);
        gD0err->SetMarkerColor(COLOR[icent]);
        gD0err->SetTitle(nameCent[icent]);
        gD0err->Write(Form("gD0_err_%s",nameCent1[icent]));
        gD0sys->SetMarkerStyle(MARKERSTYLE[0]);
        gD0sys->SetMarkerSize(2.);
        gD0sys->SetMarkerColor(COLOR[icent]);
        gD0sys->SetTitle(nameCent[icent]);
        gD0sys->Write(Form("gD0_sys_%s",nameCent1[icent]));
    }
    fout->Close();
    
}
