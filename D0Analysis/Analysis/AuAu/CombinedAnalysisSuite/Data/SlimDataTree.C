Double_t GetDeltaPhi(Double_t phia, Double_t phib)
{
    Double_t pi = TMath::Pi();

    if (phia < 0)         phia += 2*pi;
    else if (phia > 2*pi) phia -= 2*pi;
    if (phib < 0)         phib += 2*pi;
    else if (phib > 2*pi) phib -= 2*pi;
    Double_t dphi = phib - phia;

    if (dphi < -1.0*pi)  dphi += 2*pi;
    else if (dphi > pi)  dphi -= 2*pi;

    // phi is between  -pi < phi < pi                                                                                                                                                                                                                                       
    return dphi;
}

Double_t dPhi(Double_t phi1, Double_t phi2) {
  Double_t deltaPhi;
  deltaPhi = abs(phi1 - phi2); //TODO absolute values
  if (deltaPhi>(2*TMath::Pi()))  deltaPhi-=2*(TMath::Pi());
  if (deltaPhi<(0*TMath::Pi())) deltaPhi+=2*(TMath::Pi()); 

  if (deltaPhi > TMath::Pi()) deltaPhi= 2*(TMath::Pi()) - deltaPhi;
  return deltaPhi;   // dphi in [0, 2Pi]
}

void SlimDataTree(TString mode = "CS1"){
    // TFile *f_D = new TFile("/Volumes/WorkDrive/STAR-Workspace/D0Analysis/Analysis/AuAu/Files/D0Jets_WithCS_JetTree_May20.root");
    // TChain *t1 = (TChain*)f_D->Get("JetTree/D0Jets");//JetTreeSigBg/SignalPlusBackground");

  TFile *f_D;
  // if (mode == "AreaBased")f_D = new TFile("/Volumes/WorkDrive/STAR-Workspace/D0Analysis/Analysis/AuAu/Files/May27_2023/D0JetsWithAreaBasedBG.root");
  // else if (mode == "CS1") f_D = new TFile("/Volumes/WorkDrive/STAR-Workspace/D0Analysis/Analysis/AuAu/Files/May27_2023/D0JetsWithCSBasedBG.root");

  if (mode == "AreaBased")f_D = new TFile("/Volumes/WorkDrive/STAR-Workspace/D0Analysis/Analysis/AuAu/Files/Jun20_2023/D0Jets_CS_Jun20_2023.root");
  else if (mode == "CS1") f_D = new TFile("/Volumes/WorkDrive/STAR-Workspace/D0Analysis/Analysis/AuAu/Files/Jun20_2023/D0Jets_CS_Jun20_2023.root");


  TChain *t1 = (TChain*)f_D->Get("JetTree/D0Jets");//JetTreeSigBg/SignalPlusBackground");

  // TFile *f_D_Secondary = new TFile("/Volumes/WorkDrive/STAR-Workspace/D0Analysis/Analysis/AuAu/Files/May27_2023/D0JetsCombined_May28_2023.root");
  // TChain *t2 = (TChain*)f_D_Secondary->Get("JetTree/D0Jets");//JetTreeSigBg/SignalPlusBackground");

  double binning[6] = {0,0.05,0.1,0.2,0.3,0.4};

  TH1F *hD0Pt = new TH1F("hD0Pt","hD0Pt",50,0,50);
  TH1F *hJetPt = new TH1F("hJetPt","hJetPt",50,0,50);
  TH1F *hZ = new TH1F("hZ","hZ",10,0,1);
  TH1F *hR = new TH1F("hR","hR",5,binning);

  TH1F *hD0Pt_SB = new TH1F("hD0Pt_SB","hD0Pt_SB",50,0,50);
  TH1F *hJetPt_SB = new TH1F("hJetPt_SB","hJetPt_SB",50,0,50);
  TH1F *hZ_SB = new TH1F("hZ_SB","hZ_SB",10,0,1);
  TH1F *hR_SB = new TH1F("hR_SB","hR_SB",5,binning);

  TH1F *hD0Pt_Low = new TH1F("hD0Pt_Low","hD0Pt_Low",50,0,50);
  TH1F *hJetPt_Low = new TH1F("hJetPt_Low","hJetPt_Low",50,0,50);
  TH1F *hZ_Low = new TH1F("hZ_Low","hZ_Low",10,0,1);
  TH1F *hR_Low = new TH1F("hR_Low","hR_Low",5,binning);

  TH1F *hD0Pt_SB_Low = new TH1F("hD0Pt_SB_Low","hD0Pt_SB_Low",50,0,50);
  TH1F *hJetPt_SB_Low = new TH1F("hJetPt_SB_Low","hJetPt_SB_Low",50,0,50);
  TH1F *hZ_SB_Low = new TH1F("hZ_SB_Low","hZ_SB_Low",10,0,1);
  TH1F *hR_SB_Low = new TH1F("hR_SB_Low","hR_SB_Low",5,binning);

  float JetArea;
  float D0Mass;float D0Pt;float D0Phi; float D0Eta;
  float KaonPt;float KaonPhi; float KaonEta;float PionPt;float PionPhi; float PionEta;
  float JetCorrPt; float JetCSEta; float JetEta; float JetCSPhi; float JetPhi;float JetPt;
  float JetHighestTrackPt;
  float Centrality;
  float Weight;
  float TrackPt[10000];
  float TrackEta[10000];
  float TrackPhi[10000];
  float TrackID[10000];
  float TrackPx[10000];
  float TrackPy[10000];
  float TrackPz[10000];
  int JetNConst;
  float _mM; float _mPt; float _mEta; float _mJPt; float _mJPtArea; float _mR; float _mZ;float _mHPt;float _mJEta; float _mZ_Area; float _m_gRefMultCorr;
  int _mCen;float _mWeight;float _mJetArea;int RunID;int EventId; float jetptfromarea; int jetnconstfromarea;
  int _mConstArea; int _mConstCS;

  t1->SetBranchStatus("*",0);
  t1->SetBranchStatus("D0Mass",1);
  t1->SetBranchStatus("D0Pt",1);
  t1->SetBranchStatus("D0Phi",1);
  t1->SetBranchStatus("D0Eta",1);
  t1->SetBranchStatus("KaonPt",1);
  t1->SetBranchStatus("KaonPhi",1);
  t1->SetBranchStatus("KaonEta",1);
  t1->SetBranchStatus("PionPt",1);
  t1->SetBranchStatus("PionPhi",1);
  t1->SetBranchStatus("PionEta",1);
  t1->SetBranchStatus("Centrality",1);
  t1->SetBranchStatus("gRefMult", 1);
  t1->SetBranchStatus("Weight",1);
  t1->SetBranchStatus("RunID",1);
  t1->SetBranchStatus("EventId",1);

  if (mode == "AreaBased"){
    t1->SetBranchStatus("JetPt",1);
    t1->SetBranchStatus("JetEta",1);
    t1->SetBranchStatus("JetPhi",1);
    t1->SetBranchStatus("JetArea",1);
    t1->SetBranchStatus("JetNConst",1);
  }

  else if (mode == "CS1"){
    t1->SetBranchStatus("JetPt",1);
    t1->SetBranchStatus("JetCSCorrPt",1);
    t1->SetBranchStatus("JetEta",1);
    t1->SetBranchStatus("JetPhi",1);
    t1->SetBranchStatus("JetArea",1);
    t1->SetBranchStatus("JetNConst",1);
    t1->SetBranchStatus("JetCSNConst",1);
  }

  else if (mode == "CS2"){
    t1->SetBranchStatus("JetCS2CorrPt",1);
    t1->SetBranchStatus("JetEta",1);
    t1->SetBranchStatus("JetPhi",1);
    t1->SetBranchStatus("JetArea",1);
    t1->SetBranchStatus("JetCS2NConst",1);
  }

  t1->SetBranchAddress( "D0Mass" , &D0Mass );
  t1->SetBranchAddress( "D0Pt" , &D0Pt );
  t1->SetBranchAddress( "D0Eta" , &D0Eta );
  t1->SetBranchAddress( "D0Phi" , &D0Phi );
  t1->SetBranchAddress( "KaonPt" , &KaonPt );
  t1->SetBranchAddress( "KaonPhi" , &KaonPhi );
  t1->SetBranchAddress( "KaonEta" , &KaonEta );
  t1->SetBranchAddress( "PionPt" , &PionPt );
  t1->SetBranchAddress( "PionPhi" , &PionPhi );
  t1->SetBranchAddress( "PionEta" , &PionEta );
  t1->SetBranchAddress( "Centrality" , &Centrality );
  t1->SetBranchAddress( "gRefMult" , &_m_gRefMultCorr );
  t1->SetBranchAddress( "Weight" , &Weight );
  t1->SetBranchAddress( "RunID" , &RunID );
  t1->SetBranchAddress( "EventId" , &EventId );

  if (mode == "AreaBased"){
    t1->SetBranchAddress( "JetPt" , &JetCorrPt );
    t1->SetBranchAddress( "JetEta" , &JetEta );
    t1->SetBranchAddress( "JetPhi" , &JetPhi );
    t1->SetBranchAddress( "JetArea" , &JetArea );
    t1->SetBranchAddress( "JetNConst" , &JetNConst );
  }

  else if (mode == "CS1"){
    t1->SetBranchAddress( "JetPt" , &jetptfromarea );
    t1->SetBranchAddress( "JetCSCorrPt" , &JetCorrPt );
    t1->SetBranchAddress( "JetEta" , &JetEta );
    t1->SetBranchAddress( "JetPhi" , &JetPhi );
    t1->SetBranchAddress( "JetArea" , &JetArea );
    t1->SetBranchAddress( "JetNConst" , &jetnconstfromarea );
    t1->SetBranchAddress( "JetCSNConst" , &JetNConst );
  }

  else if (mode == "CS2"){
    t1->SetBranchAddress( "JetCS2CorrPt" , &JetCorrPt );
    t1->SetBranchAddress( "JetEta" , &JetEta );
    t1->SetBranchAddress( "JetPhi" , &JetPhi );
    t1->SetBranchAddress( "JetArea" , &JetArea );
    t1->SetBranchAddress( "JetNConst" , &JetNConst );
  }
  
  TFile *f_D1 = new TFile(Form("out_final_%s_Jun20.root", mode.Data()),"RECREATE");
  TTree *t = new TTree("t","t");

  t->Branch("mM",&_mM,"mM/F");
  t->Branch("mR",&_mR,"mR/F");
  t->Branch("mZ",&_mZ,"mZ/F");
  t->Branch("mZArea",&_mZ_Area,"mZArea/F");
  t->Branch("mPt",&_mPt,"mPt/F");
  t->Branch("mEta",&_mEta,"mEta/F");
  t->Branch("mJPt",&_mJPt,"mJPt/F");
  t->Branch("mJPtArea",&_mJPtArea,"mJPtArea/F");
  t->Branch("mCen",&_mCen,"mCen/I");
  t->Branch("mgRefMultCorr",&_m_gRefMultCorr,"mgRefMultCorr/F");
//   t->Branch("mHPt",&_mHPt,"mHPt/F");
  t->Branch("mJEta",&_mJEta,"mJEta/F");
  t->Branch("mWeight",&_mWeight,"mWeight/F");
  t->Branch("mJetArea",&_mJetArea,"mJetArea/F");
  t->Branch("mConstArea",&_mConstArea,"mConstArea/I");
  t->Branch("mConstCS",&_mConstCS,"mConstCS/I");

  int loop = t1->GetEntries();//   

  // int loop = 100;     
  ofstream output1 ("IDLookUp.txt",ios::app);      
  for(int i =0;i<loop;i++){
      if(i%10000==0)cout << "On "<< i << " out of " << loop << " " << float(i)/loop*100 << "%" << endl;
      t1->GetEntry(i);
      // t2->GetEntry(i);
      if(D0Mass<1.75 || D0Mass>2.02)continue;
      if(PionPt<0.6)continue;
      if(KaonPt<0.6)continue;
      if (JetNConst == 0) {
        if (mode == "CS1" && jetnconstfromarea != 0) cout << "Entry = " << i << endl;
        continue;
      }

    //   if (JetArea == 0 && JetCorrPt == 0 && JetEta == 0 && JetPhi == 0) continue;

      TVector3 p1( KaonPt * cos(KaonPhi),KaonPt * sin(KaonPhi),KaonPt * sinh(KaonEta));
      TVector3 p2( PionPt * cos(PionPhi),PionPt * sin(PionPhi),PionPt * sinh(PionEta));
      // TVector3 p3( D0Pt * cos(D0Phi), D0Pt * sin(D0Phi), D0Pt * sinh(D0Eta));
      TVector3 p3;

      p3 = p1 + p2;
      D0Pt = p3.Pt();
      D0Eta = p3.Eta();
      D0Phi = p3.Phi();
     
      double mR = sqrt((D0Eta-JetEta)*(D0Eta-JetEta)+dPhi(D0Phi,JetPhi)*dPhi(D0Phi,JetPhi));

      if (mR > 0.8){
        cout << "Event Id: " << EventId << " Run Id: " << RunID << endl;
        cout << "D0Pt: " << D0Pt << " D0Eta: " << D0Eta << " D0Phi: " << D0Phi << " DeltaR: " << mR << endl;
        cout << "JetPt: " << JetCorrPt << " JetEta: " << JetEta << " JetPhi: " << JetPhi << " JetArea: " << JetArea << endl;
      }

       if (mode == "CS1"){
        if (JetCorrPt < D0Pt) JetCorrPt = D0Pt;
       }

//      TVector3 j1( JetPt*cos(JetPhi),JetPt*sin(JetPhi),JetPt*sinh(JetEta));
      TVector3 j1( JetCorrPt*cos(JetPhi),JetCorrPt*sin(JetPhi),JetCorrPt*sinh(JetEta));    
      //double mZ = TMath::Cos(p3.Angle(j1))*p3.Mag()/j1.Mag();
      double mZ = (j1.Px()*p3.Px() + j1.Py()*p3.Py())/TMath::Power(JetCorrPt,2);
      double mZArea;
      
      int ccen = -1;
      if(Centrality == 0)ccen = 70;
      if(Centrality == 1)ccen = 60;
      if(Centrality == 2)ccen = 50;
      if(Centrality == 3)ccen = 40;
      if(Centrality == 4)ccen = 30;
      if(Centrality == 5)ccen = 20;
      if(Centrality == 6)ccen = 10;
      if(Centrality == 7)ccen = 5;
      if(Centrality == 8)ccen = 0;
      if(fabs(JetEta)>=0.6)continue;

      if (mode == "CS1"){
        // if (JetNConst == 1 && JetCorrPt < D0Pt) JetCorrPt = D0Pt;
        mZArea = p3.Pt()*TMath::Cos(p3.Phi()-JetPhi)/jetptfromarea;
      }

      // if (mode == "CS1"){
      //   if (jetptfromarea < 3.) continue; // This is to bring parity with the area based bg measurement
      // }
      _mCen = ccen;//Centrality;      
      _mM = D0Mass;
      _mR = mR;
      _mZ = mZ;
      _mZ_Area = mZArea;
      _mJEta =JetEta;
      _mPt = D0Pt;
      _mEta = p3.PseudoRapidity();
      _mJPt = JetCorrPt;
      _mJPtArea = jetptfromarea;
      _mConstArea = jetnconstfromarea;
      _mConstCS = JetNConst;
    //   _mHPt = JetHighestTrackPt;
      _mWeight = Weight;
      _mJetArea = JetArea;
      t->Fill();
  }

  t->Write();
}
