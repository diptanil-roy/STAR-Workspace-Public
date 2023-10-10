void ReadTree(TString filename = "PYTHIA_BBBar_FullRange.root", TString outfilename = "BBBar.root"){
  TFile *f = new TFile(filename);
  TTree *t = (TTree *)f->Get("Jets");
  TFile *fout = new TFile(outfilename, "RECREATE");

  const int njpt_gen_bins_var = 8; //x2 MC
  double jetpt_var_bin[njpt_gen_bins_var+1] = {1,3,5,7,9,11,13,15,20};

  const int nz_gen_bins = 10; //x2 MC
  double z_gen_bin[nz_gen_bins+1] = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.};
  
  TH2D *ExternalTruthHistogram = new TH2D("ExternalTruthHistogram", "ExternalTruthHistogram", njpt_gen_bins_var, jetpt_var_bin, nz_gen_bins, z_gen_bin);

  float jetpt;
  float jeteta;
  float jetphi;
  float D0pt;
  float D0eta;
  float D0phi;
  float z;
  t->SetBranchAddress("JetPt", &jetpt);
  t->SetBranchAddress("JetEta", &jeteta);
  t->SetBranchAddress("JetPhi", &jetphi);
  t->SetBranchAddress("D0Pt", &D0pt);
  t->SetBranchAddress("D0Eta", &D0eta);
  t->SetBranchAddress("D0Phi", &D0phi);
  t->SetBranchAddress("D0Z", &z);
  for (int i = 0; i < t->GetEntries(); i++){
    t->GetEntry(i);
    if (D0pt < 1.0 || D0pt > 10.0) continue;
    if (jetpt > 20.) continue;
    ExternalTruthHistogram->Fill(jetpt, z);
  }
  fout->cd();
  ExternalTruthHistogram->Write();
  fout->Close();
}