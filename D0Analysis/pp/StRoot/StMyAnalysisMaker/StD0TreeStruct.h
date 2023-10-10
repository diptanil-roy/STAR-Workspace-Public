#pragma link C++ class vector<float> +;
#pragma link C++ class vector<vector<float> >+;
#pragma link C++ class vector<int> +;
#pragma link C++ class vector<vector<int> >+;

#ifndef StD0TreeStruct_h
#define StD0TreeStruct_h

struct StD0TreeStruct
{

    float Vz;
    float VzVPD;
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
};

#endif