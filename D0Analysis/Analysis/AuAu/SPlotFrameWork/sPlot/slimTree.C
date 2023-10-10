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
void slimTree(){
    //TFile *f_D = new TFile("/gpfs01/star/pwg_tasks/jetcorr03/Mar24_AllUpdatedFiles/DataWithWiderEtaCut.root");
    TFile *f_D = new TFile("../../Files/NewDataSample.root");
    TChain *t1 = (TChain*)f_D->Get("JetTree_Standard_UnlikeSign/D0Jets");//JetTreeSigBg/SignalPlusBackground");
    TChain *t2 = (TChain*)f_D->Get("JetTree_Standard_LikeSign/D0Jets");//JetTreeULBg/UnlikeBackground");

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
  float D0Mass;float KaonPt;float KaonPhi; float KaonEta;float PionPt;float PionPhi; float PionEta;
  float JetCorrPt; float JetEta; float JetPhi;float JetPt;
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
  float _mM; float _mPt; float _mEta; float _mJPt; float _mR; float _mZ;float _mHPt;float _mJEta; 
  int _mCen;float _mWeight;float _mJetArea;
  t1->SetBranchAddress( "D0Mass" , &D0Mass );
  t1->SetBranchAddress( "KaonPt" , &KaonPt );
  t1->SetBranchAddress( "KaonPhi" , &KaonPhi );
  t1->SetBranchAddress( "KaonEta" , &KaonEta );
  t1->SetBranchAddress( "PionPt" , &PionPt );
  t1->SetBranchAddress( "PionPhi" , &PionPhi );
  t1->SetBranchAddress( "PionEta" , &PionEta );
  t1->SetBranchAddress( "JetCorrPt" , &JetCorrPt );
  t1->SetBranchAddress( "JetPt" , &JetPt );
  t1->SetBranchAddress( "JetEta" , &JetEta );
  t1->SetBranchAddress( "JetPhi" , &JetPhi );
  t1->SetBranchAddress( "Centrality" , &Centrality );
  t1->SetBranchAddress( "TrackPx" , TrackPx );
  t1->SetBranchAddress( "TrackPy" , TrackPy );
  t1->SetBranchAddress( "TrackPz" , TrackPz );
  t1->SetBranchAddress( "TrackID" , TrackID );
  t1->SetBranchAddress( "JetNConst" , &JetNConst );
  t1->SetBranchAddress( "JetHighestTrackPt" , &JetHighestTrackPt );
  t1->SetBranchAddress( "Weight" , &Weight );
  t1->SetBranchAddress( "JetArea" , &JetArea );
  int RunID;
  int EventId;
  t1->SetBranchAddress( "RunID" , &RunID );
  t1->SetBranchAddress( "EventId" , &EventId );
  
  TFile *f_D1 = new TFile("out_final.root","RECREATE");
  TTree *t = new TTree("t","t");

  t->Branch("mM",&_mM,"mM/F");
  t->Branch("mR",&_mR,"mR/F");
  t->Branch("mZ",&_mZ,"mZ/F");
  t->Branch("mPt",&_mPt,"mPt/F");
  t->Branch("mEta",&_mEta,"mEta/F");
  t->Branch("mJPt",&_mJPt,"mJPt/F");
  t->Branch("mCen",&_mCen,"mCen/I");
  t->Branch("mHPt",&_mHPt,"mHPt/F");
  t->Branch("mJEta",&_mJEta,"mJEta/F");
  t->Branch("mWeight",&_mWeight,"mWeight/F");
  t->Branch("mJetArea",&_mJetArea,"mJetArea/F");
  int loop = t1->GetEntries();//        
  ofstream output1 ("IDLookUp.txt",ios::app);      
  for(int i =0;i<loop;i++){
      if(i%10000==0)cout << "On "<< i << " out of " << loop << " " << float(i)/loop*100 << "%" << endl;
      t1->GetEntry(i);
      if(D0Mass<1.75 || D0Mass>2.02)continue;
      if(PionPt<0.6)continue;
      if(KaonPt<0.6)continue;
      TVector3 p1( KaonPt * cos(KaonPhi),KaonPt * sin(KaonPhi),KaonPt * sinh(KaonEta));
      TVector3 p2( PionPt * cos(PionPhi),PionPt * sin(PionPhi),PionPt * sinh(PionEta));
      TVector3 p3;// = p1+p2;
      for(int idx = 0; idx<JetNConst;idx++){
	  //if(TrackID[idx]==99999){
	  if(TrackID[idx]==421){
	      p3.SetXYZ(TrackPx[idx],TrackPy[idx],TrackPz[idx]);
	      break;
	  }
      }
      double mR = sqrt((p3.PseudoRapidity()-JetEta)*(p3.PseudoRapidity()-JetEta)+GetDeltaPhi(p3.Phi(),JetPhi)*GetDeltaPhi(p3.Phi(),JetPhi));
//      TVector3 j1( JetPt*cos(JetPhi),JetPt*sin(JetPhi),JetPt*sinh(JetEta));
      TVector3 j1( JetCorrPt*cos(JetPhi),JetCorrPt*sin(JetPhi),JetCorrPt*sinh(JetEta));    
      //double mZ = TMath::Cos(p3.Angle(j1))*p3.Mag()/j1.Mag();
      double mZ = (j1.Px()*p3.Px() + j1.Py()*p3.Py())/TMath::Power(JetCorrPt,2);
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
      _mCen = ccen;//Centrality;
      _mM = D0Mass;
      _mR = mR;
      _mZ = mZ;
      _mJEta =JetEta;
      _mPt = p3.Pt();
      _mEta = p3.PseudoRapidity();
      _mJPt = JetCorrPt;
      _mHPt = JetHighestTrackPt;
      _mWeight = Weight;
      _mJetArea = JetArea;
      t->Fill();
      if(p3.Pt() > 5. )output1 << RunID << " " << EventId << endl;
      if(JetCorrPt>5.){
    	  hD0Pt->Fill(p3.Pt());
    	  hJetPt->Fill(JetCorrPt);
    	  hZ->Fill(mZ);
    	  hR->Fill(mR);
    	  if(p3.Pt()<3.){
    	      hD0Pt_Low->Fill(p3.Pt());
    	      hJetPt_Low->Fill(JetCorrPt);
    	      hZ_Low->Fill(mZ);
    	      hR_Low->Fill(mR);
    	  }

      }

  }
/*  cout <<"> Starting Second Loop" << endl;
  t2->SetBranchAddress( "D0Mass" , &D0Mass );
  t2->SetBranchAddress( "KaonPt" , &KaonPt );
  t2->SetBranchAddress( "KaonPhi" , &KaonPhi );
  t2->SetBranchAddress( "KaonEta" , &KaonEta );
  t2->SetBranchAddress( "PionPt" , &PionPt );
  t2->SetBranchAddress( "PionPhi" , &PionPhi );
  t2->SetBranchAddress( "PionEta" , &PionEta );
  t2->SetBranchAddress( "JetCorrPt" , &JetCorrPt );
  t2->SetBranchAddress( "JetPt" , &JetPt );
  t2->SetBranchAddress( "JetEta" , &JetEta );
  t2->SetBranchAddress( "JetPhi" , &JetPhi );
  t2->SetBranchAddress( "Centrality" , &Centrality );
  t2->SetBranchAddress( "TrackPx" , TrackPx );
  t2->SetBranchAddress( "TrackPy" , TrackPy );
  t2->SetBranchAddress( "TrackPz" , TrackPz );
  t2->SetBranchAddress( "TrackID" , TrackID );
  t2->SetBranchAddress( "JetNConst" , &JetNConst );
  t2->SetBranchAddress( "JetHighestTrackPt" , &JetHighestTrackPt );
  loop = t2->GetEntries();//                                                                                                                                                                                                                                         
  for(int i =0;i<loop;i++){
      if(i%10000==0)cout << "On "<< i << " out of " << loop << " " << float(i)/loop*100 << "%" << endl;
      t2->GetEntry(i);
      if(D0Mass<1.75 || D0Mass>2.02)continue;
      TVector3 p1( KaonPt * cos(KaonPhi),KaonPt * sin(KaonPhi),KaonPt * sinh(KaonEta));
      TVector3 p2( PionPt * cos(PionPhi),PionPt * sin(PionPhi),PionPt * sinh(PionEta));
      TVector3 p3;// = p1+p2;
      for(int idx = 0; idx<JetNConst;idx++){
          if(TrackID[idx]==99999){
	      p3.SetXYZ(TrackPx[idx],TrackPy[idx],TrackPz[idx]);
	      break;
	  }
      }
      double mR = sqrt((p3.PseudoRapidity()-JetEta)*(p3.PseudoRapidity()-JetEta)+GetDeltaPhi(p3.Phi(),JetPhi)*GetDeltaPhi(p3.Phi(),JetPhi));

      TVector3 j1( JetCorrPt*cos(JetPhi),JetCorrPt*sin(JetPhi),JetCorrPt*sinh(JetEta));
      double mZ = TMath::Cos(p3.Angle(j1))*p3.Mag()/j1.Mag();
      _mCen = Centrality;
      _mM = D0Mass;
      _mR = mR;
      _mZ = mZ;
      _mJEta =JetEta;
      _mPt = p3.Pt();
      _mEta = p3.PseudoRapidity();
      _mJPt = JetCorrPt;
      _mHPt = JetHighestTrackPt;
      t->Fill();
      if(JetCorrPt>5.){
	  hD0Pt_SB->Fill(p3.Pt());
	  hJetPt_SB->Fill(JetCorrPt);
	  hZ_SB->Fill(mZ);
	  hR_SB->Fill(mR);
	  if(p3.Pt()<3.){
              hD0Pt_SB_Low->Fill(p3.Pt());
              hJetPt_SB_Low->Fill(JetCorrPt);
              hZ_SB_Low->Fill(mZ);
              hR_SB_Low->Fill(mR);
          }

      }
  }
*/

  t->Write();
  hD0Pt_SB->Write();
  hJetPt_SB->Write();
  hZ_SB->Write();
  hR_SB->Write();
  hD0Pt->Write();
  hJetPt->Write();
  hZ->Write();
  hR->Write();

  hD0Pt_SB_Low->Write();
  hJetPt_SB_Low->Write();
  hZ_SB_Low->Write();
  hR_SB_Low->Write();
  hD0Pt_Low->Write();
  hJetPt_Low->Write();
  hZ_Low->Write();
  hR_Low->Write();
}
