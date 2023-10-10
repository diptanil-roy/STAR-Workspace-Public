#ifndef MuHijing_h
#define MuHijing_h

// const int JetMax = 1000;
const int ConstMax = 5000;
struct MuHijing
{
    
    // Event Data Info
    // int eventid;
    int numberofbinary;
    float impactparameter;

    int nparticles;
    float mPt[ConstMax];
    float mEta[ConstMax];
    float mPhi[ConstMax];
    float mEnergy[ConstMax];
    int mPID[ConstMax];
};

#endif