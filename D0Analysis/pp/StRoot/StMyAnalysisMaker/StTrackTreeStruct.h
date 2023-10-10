#pragma link C++ class vector<float> +;
#pragma link C++ class vector<vector<float> >+;
#pragma link C++ class vector<int> +;
#pragma link C++ class vector<vector<int> >+;

#ifndef StTrackTreeStruct_h
#define StTrackTreeStruct_h

struct StTrackTreeStruct
{
    int eventid;
    vector <double> primaryvertex;
    float dca;
    float pt;
    float eta;
    float phi;
    int charge;

    int nhitsfit;
    int nhitsmax;
    float nsigmapion;
    float nsigmakaon;
    float nsigmaproton;

    float ntofbeta;
    float ntofylocal;
    float isbemctrack;
};

#endif