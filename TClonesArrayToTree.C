class Track: public TObject{
public:
    int fIndex;
    float fPx;
    float fPy;
    float fPz;
    float fE;
    float fCharge;
    int fPID;
    Track() : fIndex(0), fPx(0), fPy(0), fPz(0), fE(0), fCharge(0), fPID(0) {}
ClassDef(Track, 1)
};


class Event : public TObject{
public:
    Event(){}
    Int_t fNtrack = 0; // number of tracks in this event
    vector <Track> fTracks; //->array of tracks
    // fTracks.clear();
    virtual ~Event() { fTracks.clear(); }
    ClassDef(Event, 1)
};

ClassImp(Track)
ClassImp(Event)

void TClonesArrayToTree(){
    TFile *f = new TFile("TClone.root", "RECREATE");
    TTree *tree = new TTree("T", "Testing TClonesArrays");
    
    TClonesArray *a = new TClonesArray("Event");
    TClonesArray &aa = *a;

    Event *e = new Event();
    Event &ee = *e;
    // tree->Branch("test", &a);

    tree->Branch("Events", "Event", &e);

    // create 1000 TLorentzVectors in a loop.
    for (Int_t i=0;i<1000;i++) {
        
        ee.fNtrack = 100;
        ee.fTracks.clear();
        // create 10 tracks in a loop

        // cout << ee.fNtrack << endl;
        
        for (Int_t j=0;j<ee.fNtrack;j++) {
            Track *t = new Track();
            t->fIndex = j;
            t->fPx = gRandom->Uniform(0,1);
            t->fPy = gRandom->Uniform(0,1);
            t->fPz = gRandom->Uniform(0,1);
            t->fE = gRandom->Uniform(0,1);
            t->fCharge = gRandom->Uniform(0,1);
            t->fPID = j;

            ee.fTracks.push_back(*t);

            // cout << "Track " << j << " px: " << t->fPx << endl;
        }

        // new(aa[i]) Event(e);

        tree->Fill();

        // delete e;
    }

    // fill tree... (likely to be in a loop)
    

    // write tree to file...
    f->Write();
    f->Close();
    //delete f;
}