#pragma link C++ class vector<float> +;
#pragma link C++ class vector<vector<float> >+;
#pragma link C++ class vector<int> +;
#pragma link C++ class vector<vector<int> >+;

#ifndef StFakeJetTree_h
#define StFakeJetTree_h

// const int JetMax = 1000;
const int MaxTrack = 5000;
struct StFakeJetTree
{
    /* There's only one jet per event, due to the structure of our analysis. If an event has more than one jet of interest, that eventid will have a frequency > 1.
    Nothing else should change. 
    The event level info we store are            : run id, event id, event refmult, event centrality, event triggers.
    The jet level info we store are              : jetpt, jetcorrectedpt, jeteta, jetphi, jetarea, jetradius, jetenergy, jetnef (neutral energy fraction), fRhoVal used, 
                                                   highest energy track pt, number of constituents
    The constituent level info we store are      : trackid, trackpt, tracketa, trackphi, trackpx, trackpy, trackpz                                 


    */
    // Jet Data Info
    int runid;
    int eventid;
    float refmult;
    float centrality;
    vector <unsigned int> triggers;
    vector <double> primaryvertex;
    vector <double> primaryvertexerror;

    float hardpt;
    float hardeta;
    float hardphi;

    float jetpt;
    float jetcorrectedpt;
    float jeteta;
    float jetphi;
    float jetarea;
    float jetradius;
    float jetenergy;
    float jetnef;
    float fRhoValforjet;
    float jethighesttrackpt;
    int numberofconstituents;

    float mTrackID[MaxTrack];
    float mTrackPt[MaxTrack];
    float mTrackEta[MaxTrack];
    float mTrackPhi[MaxTrack];
    float mTrackPx[MaxTrack];
    float mTrackPy[MaxTrack];
    float mTrackPz[MaxTrack];
    float mTrackCharge[MaxTrack];

};

#endif