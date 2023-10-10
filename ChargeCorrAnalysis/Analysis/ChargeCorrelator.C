TH1D *ReturnRC(TH1 *LS, TH1* US){
    TH1::SetDefaultSumw2();

    TH1D *tmp1 = (TH1D*)LS->Clone();
    TH1D *tmp2 = (TH1D*)US->Clone();
    TH1D *tmp3 = (TH1D*)LS->Clone();

    tmp3->Add(US);
    tmp1->Add(tmp2, -1);

    tmp1->Divide(tmp3);

    TString name = LS->GetName();
    name.ReplaceAll("LS", "RC");

    tmp1->SetNameTitle(name.Data(), name.Data());
    return tmp1;
}

void ChargeCorrelator(){

    TH1::SetDefaultSumw2();

    TFile *f = new TFile("../Files/JetTree_ChargeCorr.root", "READ");
    f->cd("JetTree");

    TTree *t = (TTree*)gDirectory->Get("ChargeCorrelation");
    cout << t->GetName() << endl;

    int JetNConst;
    float Centrality, JetPt, LeadTrackPt, SubLeadTrackPt, LeadTrackPhi, SubLeadTrackPhi, LeadTrackEta, SubLeadTrackEta;
    bool LeadTrackIsPion, LeadTrackIsKaon;
    bool SubLeadTrackIsPion, SubLeadTrackIsKaon;
    int LeadTrackCharge, SubLeadTrackCharge;

    t->SetBranchStatus("*", 0);
    t->SetBranchStatus("Centrality", 1);
    t->SetBranchStatus("JetPt", 1);
    t->SetBranchStatus("JetNConst", 1);
    t->SetBranchStatus("LeadTrackPt", 1);
    t->SetBranchStatus("LeadTrackEta", 1);
    t->SetBranchStatus("LeadTrackPhi", 1);
    t->SetBranchStatus("SubLeadTrackPt", 1);
    t->SetBranchStatus("SubLeadTrackEta", 1);
    t->SetBranchStatus("SubLeadTrackPhi", 1);
    t->SetBranchStatus("LeadTrackIsPion", 1);
    t->SetBranchStatus("LeadTrackIsKaon", 1);
    t->SetBranchStatus("SubLeadTrackIsPion", 1);
    t->SetBranchStatus("SubLeadTrackIsKaon", 1);
    t->SetBranchStatus("LeadTrackCharge", 1);
    t->SetBranchStatus("SubLeadTrackCharge", 1);

    t->SetBranchAddress("Centrality", &Centrality);
    t->SetBranchAddress("JetPt", &JetPt);
    t->SetBranchAddress("JetNConst", &JetNConst);
    t->SetBranchAddress("LeadTrackPt", &LeadTrackPt);
    t->SetBranchAddress("LeadTrackEta", &LeadTrackEta);
    t->SetBranchAddress("LeadTrackPhi", &LeadTrackPhi);
    t->SetBranchAddress("SubLeadTrackPt", &SubLeadTrackPt);
    t->SetBranchAddress("SubLeadTrackEta", &SubLeadTrackEta);
    t->SetBranchAddress("SubLeadTrackPhi", &SubLeadTrackPhi);
    t->SetBranchAddress("LeadTrackIsPion", &LeadTrackIsPion);
    t->SetBranchAddress("LeadTrackIsKaon", &LeadTrackIsKaon);
    t->SetBranchAddress("SubLeadTrackIsPion", &SubLeadTrackIsPion);
    t->SetBranchAddress("SubLeadTrackIsKaon", &SubLeadTrackIsKaon);
    t->SetBranchAddress("LeadTrackCharge", &LeadTrackCharge);
    t->SetBranchAddress("SubLeadTrackCharge", &SubLeadTrackCharge);

    // Define the histograms above for 3 centrality classes. Use arrays of histogram pointers
    // to make the code more compact.

    const int numCentralityClasses = 4; // 0 - 80%, 0 - 10%, 10 - 50%, 50 - 80%
    TH1D* LSVsJetMultiplicity[numCentralityClasses];
    TH1D* USVsJetMultiplicity[numCentralityClasses];
    TH1D* LSVsJetPt[numCentralityClasses];
    TH1D* USVsJetPt[numCentralityClasses];
    TH1D* LSVsZ[numCentralityClasses];
    TH1D* USVsZ[numCentralityClasses];
    TH1D* LSVsKPerp[numCentralityClasses];
    TH1D* USVsKPerp[numCentralityClasses];

    for (int i = 0; i < numCentralityClasses; i++) {
        // Jet Multiplicity histograms
        LSVsJetMultiplicity[i] = new TH1D(Form("LSVsJetMultiplicity_%d", i), Form("LSVsJetMultiplicity_%d", i), 10, 0.5, 10.5);
        USVsJetMultiplicity[i] = new TH1D(Form("USVsJetMultiplicity_%d", i), Form("USVsJetMultiplicity_%d", i), 10, 0.5, 10.5);
    
        // Jet Pt histograms
        LSVsJetPt[i] = new TH1D(Form("LSVsJetPt_%d", i), Form("LSVsJetPt_%d", i), 20, 10, 50);
        USVsJetPt[i] = new TH1D(Form("USVsJetPt_%d", i), Form("USVsJetPt_%d", i), 20, 10, 50);
    
        // Z histograms
        LSVsZ[i] = new TH1D(Form("LSVsZ_%d", i), Form("LSVsZ_%d", i), 20, 0, 1);
        USVsZ[i] = new TH1D(Form("USVsZ_%d", i), Form("USVsZ_%d", i), 20, 0, 1);
    
        // KPerp histograms
        LSVsKPerp[i] = new TH1D(Form("LSVsKPerp_%d", i), Form("LSVsKPerp_%d", i), 20, 0, 5);
        USVsKPerp[i] = new TH1D(Form("USVsKPerp_%d", i), Form("USVsKPerp_%d", i), 20, 0, 5);
    }

    int nentries = t->GetEntries();

    cout << "Total Entries: " << nentries << endl;

    for (int i = 0; i < nentries; i++){
        t->GetEntry(i);

        if (JetNConst < 2) continue;

        if (JetPt < 10) continue;

        int centralityClass = -99;

        if (Centrality < 10) centralityClass = 1;
        else if (Centrality < 50) centralityClass = 2;
        else if (Centrality < 80) centralityClass = 3;
        
        TVector3 LeadTrack, SubLeadTrack;
        LeadTrack.SetPtEtaPhi(LeadTrackPt, LeadTrackEta, LeadTrackPhi);
        SubLeadTrack.SetPtEtaPhi(SubLeadTrackPt, SubLeadTrackEta, SubLeadTrackPhi);

        double z = SubLeadTrack.Pt()/(LeadTrack + SubLeadTrack).Pt();
        double deltar = TMath::Sqrt(pow(LeadTrack.Eta() - SubLeadTrack.Eta(), 2) + pow(LeadTrack.Angle(SubLeadTrack), 2));
        double kperp = SubLeadTrack.Pt()*deltar;
        
        if (LeadTrackCharge == SubLeadTrackCharge){
            LSVsJetMultiplicity[0]->Fill(JetNConst);
            LSVsJetPt[0]->Fill(JetPt);
            LSVsZ[0]->Fill(z);
            LSVsKPerp[0]->Fill(kperp);

            LSVsJetMultiplicity[centralityClass]->Fill(JetNConst);
            LSVsJetPt[centralityClass]->Fill(JetPt);
            LSVsZ[centralityClass]->Fill(z);
            LSVsKPerp[centralityClass]->Fill(kperp);
        }
        else{
            USVsJetMultiplicity[0]->Fill(JetNConst);
            USVsJetPt[0]->Fill(JetPt);
            USVsZ[0]->Fill(z);
            USVsKPerp[0]->Fill(kperp);

            USVsJetMultiplicity[centralityClass]->Fill(JetNConst);
            USVsJetPt[centralityClass]->Fill(JetPt);
            USVsZ[centralityClass]->Fill(z);
            USVsKPerp[centralityClass]->Fill(kperp);
        }
    }

    TH1D *RCJetMultiplicity[numCentralityClasses];
    TH1D *RCJetPt[numCentralityClasses];
    TH1D *RCZ[numCentralityClasses];
    TH1D *RCKPerp[numCentralityClasses];

    for (int i = 0; i < numCentralityClasses; i++){
        RCJetMultiplicity[i] = ReturnRC(LSVsJetMultiplicity[i], USVsJetMultiplicity[i]);
        RCJetPt[i] = ReturnRC(LSVsJetPt[i], USVsJetPt[i]);
        RCZ[i] = ReturnRC(LSVsZ[i], USVsZ[i]);
        RCKPerp[i] = ReturnRC(LSVsKPerp[i], USVsKPerp[i]);
    }

    TFile *fout = new TFile("ChargeCorrelator_Output.root", "RECREATE");
    fout->cd();

    for (int i = 0; i < numCentralityClasses; i++){
        LSVsJetMultiplicity[i]->Write();
        USVsJetMultiplicity[i]->Write();
        RCJetMultiplicity[i]->Write();
        LSVsJetPt[i]->Write();
        USVsJetPt[i]->Write();
        RCJetPt[i]->Write();
        LSVsZ[i]->Write();
        USVsZ[i]->Write();
        RCZ[i]->Write();
        LSVsKPerp[i]->Write();
        USVsKPerp[i]->Write();
        RCKPerp[i]->Write();
    }

    fout->Close();
}