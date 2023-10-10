#pragma link C++ class vector<float> +;
#pragma link C++ class vector<vector<float> >+;
#pragma link C++ class vector<int> +;
#pragma link C++ class vector<vector<int> >+;

#ifndef StJetTreeStruct_h
#define StJetTreeStruct_h

// const int JetMax = 1000;
const int ConstMax = 5000;
struct StJetTreeStruct
{
    /* There's only one jet per event, due to the structure of our analysis. If an event has more than one jet of interest, that eventid will have a frequency > 1.
    Nothing else should change. 
    The event level info we store are            : run id, event id, event refmult, event centrality, event triggers.
    The jet level info we store are              : jetpt, jetcorrectedpt, jeteta, jetphi, jetarea, jetradius, jetenergy, jetnef (neutral energy fraction), fRhoVal used, 
                                                   highest energy track pt, number of constituents
    The constituent level info we store are      : trackid, trackpt, tracketa, trackphi, trackpx, trackpy, trackpz                                 


    */
    // Jet Data Info
    float Vz;
    float VzVPD;
    vector <unsigned int> Triggers;
    float D0mass;
    float D0pt;
    float D0eta;
    float D0phi;

    float pionpt;
    float pioneta;
    float pionphi;
    float piondca;
    float pionnhitsfit;
    float pionnsigmapion;
    float pionnsigmakaon;
    float piontofbeta;       // Save out (1/beta - 1/betapion)
    float ntofpionylocal; // Save out ylocal information

    float kaonpt;
    float kaoneta;
    float kaonphi;
    float kaondca;
    float kaonnhitsfit;
    float kaonnsigmapion;
    float kaonnsigmakaon;
    float kaontofbeta;       // Save out (1/beta - 1/betapion)
    float ntofkaonylocal; // Save out ylocal information

    float jetpt;
    float jeteta;
    float jetphi;
    float jetarea;
    float jetenergy;
    float jetnef;
    int numberofconstituents;

    float mTrackPt[ConstMax];
    float mTrackEta[ConstMax];
    float mTrackPhi[ConstMax];
    float mTrackCharge[ConstMax];

};

#endif
