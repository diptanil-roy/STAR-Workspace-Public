void test2(){
    TFile *f = new TFile("test.root","RECREATE");
    TTree *t = new TTree("t","t");

    int nparticles = 0;
    float mPt[5000];

    t->Branch("NParticles", &nparticles, "NParticles/I");
    t->Branch("TrackPt", mPt, "TrackPt[NParticles]/F");

    for (int i = 0; i < 10; i++){
        nparticles = gRandom->Integer(10);

        // Reset mPt array to 0 for each event
        for (int j = 0; j < 5000; j++){
            mPt[j] = 0.0;
        }
        
        for (int j = 0; j < nparticles; j++){
            mPt[j] = gRandom->Uniform(0,10);
        }

        t->Fill();
    }
    t->Write();
    f->Close();
}