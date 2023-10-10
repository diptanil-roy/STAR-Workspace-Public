double fillHistos(char file[100], TH2F* PT2D);
void checkBKG(){
  gROOT->ProcessLine(".x ~/myStyle.C");
  TH1::SetDefaultSumw2();

  TFile *inf = new TFile("Histograms.root","READ");
  D0Pt=(TH1F*)inf->Get("D0Pt");
  D0Pt_1=(TH1F*)inf->Get("D0Pt_0_10");
  D0Pt_2=(TH1F*)inf->Get("D0Pt_10_40");
  D0Pt_3=(TH1F*)inf->Get("D0Pt_40_80");
  
  JetPt=(TH1F*)inf->Get("JetPt");
  JetPt_1=(TH1F*)inf->Get("JetPt_0_10");
  JetPt_2=(TH1F*)inf->Get("JetPt_10_40");
  JetPt_3=(TH1F*)inf->Get("JetPt_40_80");
  
  Z=(TH1F*)inf->Get("Z");
  Z_1=(TH1F*)inf->Get("Z_0_10");
  Z_2=(TH1F*)inf->Get("Z_10_40");
  Z_3=(TH1F*)inf->Get("Z_40_80");
  
  R=(TH1F*)inf->Get("R");
  R_1=(TH1F*)inf->Get("R_0_10");
  R_2=(TH1F*)inf->Get("R_10_40");
  R_3=(TH1F*)inf->Get("R_40_80");

  double sumw=0;
  TH1F *hD0 = new TH1F("hD0","hD0",50,0,50);
  TH1F *hJet = new TH1F("hJet","hJet",50,0,50);
  sumw+=fillHistos("../Production/output/output_1/all.root",hJet,hD0,D0Pt_3);  
  

  TCanvas *c1 = new TCanvas("c1","c1");
  JetPt_3->DrawNormalized();
  hJet->DrawNormalized("hist same");

  TCanvas *c2 = new TCanvas("c2","c2");
  D0Pt_3->DrawNormalized();
  hD0->DrawNormalized("hist same");
}



double fillHistos(char file[100], TH1F* PT2D, TH1F* d0,TH1F* WW){
    int sumw=0;
    TFile *f_D = new TFile(file);
    TChain* ch = (TChain*) f_D->Get("JetTreeEmbed/Embed");
    float gpt; float rpt;float  mCen;
    ch->SetBranchAddress( "PionPt" , &gpt );
    ch->SetBranchAddress( "JetCorrPt" , &rpt );
    ch->SetBranchAddress( "Centrality" , &mCen );
    int loop = ch->GetEntries();
    for(int i =0;i<loop;i++){
        if(i%10000==0)cout << "On "<< i << " out of " << loop << " " << float(i)/loop*100 << "%" << endl;
        ch->GetEntry(i);
        sumw++;
	//if(gpt<5)continue;
        if(rpt<5)continue;
	double weight = WW->GetBinContent(WW->FindBin(gpt));
        if(mCen>=40 && mCen<80){
	    PT2D->Fill(rpt,weight);
	    d0->Fill(gpt,weight);
        }
    }
    f_D->Close();
    return sumw;
}
