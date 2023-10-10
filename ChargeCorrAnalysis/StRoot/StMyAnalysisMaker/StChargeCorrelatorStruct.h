#ifndef StChargeCorrelatorStruct_h
#define StChargeCorrelatorStruct_h

struct StChargeCorrelatorStruct
{                             
    // Jet Data Info
    int runid;
    int eventid;
    float refmult;
    float grefmult;
    float centrality;
    float weight;
    bool isMB;
    bool isHT1;
    bool isHT2;
    bool isHT3;

    float jetpt;
    float jeteta;
    float jetphi;
    float jetarea;
    float jetenergy;
    float jetnef;
    int numberofconstituents;

    float leadtrackpt;
    float leadtracketa;
    float leadtrackphi;
    int leadtrackcharge;
    bool leadtrackispion;
    bool leadtrackiskaon;
    bool leadtrackisproton;

    float subleadtrackpt;
    float subleadtracketa;
    float subleadtrackphi;
    int subleadtrackcharge;
    bool subleadtrackispion;
    bool subleadtrackiskaon;
    bool subleadtrackisproton;
};

#endif