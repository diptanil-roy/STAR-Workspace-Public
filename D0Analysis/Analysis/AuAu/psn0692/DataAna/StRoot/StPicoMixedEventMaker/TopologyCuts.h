#ifndef TopologyCuts__H
#define TopologyCuts__H
#include <vector>

//#include "TopologyCuts.h"
//#include "StMixerCuts.h"
/* **************************************************
 *  Topology Cut Sets
 *
 *  Authors:  Michael Lomnitz (mrlomnitz@lbl.gov)
 *
 * **************************************************
 */
namespace topoCuts
{
    struct TopologicalCuts
    {
        float  RapidityCut;
        int    nPtBins;
        int    nCentBins;
        std::vector<float> PtEdge;
        std::vector<float> CentEdge;
        std::vector< std::vector<float> > dcaV0ToPv;
        std::vector< std::vector<float> > decayLength;
        std::vector< std::vector<float> > cosTheta;
        std::vector< std::vector<float> > dcaDaughters;
        std::vector< std::vector<float> > kDca;
        std::vector< std::vector<float> > pDca;
        TopologicalCuts(float const RapidityCut1,
                        const int   nPtBins1,
                        const float PtEdge1[],
                        const int   nCentBins1,
                        const float CentEdge1[],
                        const float dcaV0ToPv1[][6],
                        const float decayLength1[][6],
                        //const float cosTheta1[][6],
                        const float dcaDaughters1[][6],
                        const float kDca1[][6],
                        const float pDca1[][6])
                        :
        RapidityCut(RapidityCut1),
        nPtBins(nPtBins1),
        nCentBins(nCentBins1)
        {
            for(int ipt=0; ipt<nPtBins+1; ipt++) PtEdge.push_back(PtEdge1[ipt]);
            for(int icent=0; icent<nCentBins+1; icent++) CentEdge.push_back(CentEdge1[icent]);
            for(int icent=0; icent<nCentBins1; icent++) {
                std::vector<float> tmpDcaK, tmpDcaP, tmpDcaD0, tmpDca12, tmpDecayL, tmpCosTheta;
                for(int ipt=0; ipt<nPtBins1; ipt++) {
                    tmpDcaD0.push_back(dcaV0ToPv1[icent][ipt]);
                    tmpDecayL.push_back(decayLength1[icent][ipt]);
                    tmpCosTheta.push_back(0.95);
                    tmpDca12.push_back(dcaDaughters1[icent][ipt]);
                    tmpDcaK.push_back(kDca1[icent][ipt]);
                    tmpDcaP.push_back(pDca1[icent][ipt]);
                }
                dcaV0ToPv.push_back(tmpDcaD0);
                decayLength.push_back(tmpDecayL);
                dcaDaughters.push_back(tmpDca12);
                kDca.push_back(tmpDcaK);
                pDca.push_back(tmpDcaP);
                cosTheta.push_back(tmpCosTheta);
            }
         };
    };
    
    const int nCent = 5;
    const float CentEdge[nCent+1] = { -0.5, 1.5, 3.5, 5.5, 6.5, 8.5 };
    const char nameCent[nCent][100] = {"60-80%", "40-60%", "20-40%", "10-20%", "0-10%"};
    
    const int nPtBins = 6;
    const float PtEdge[nPtBins+1] = {0, 0.5, 1., 2., 3., 5., 10.};
    
    const float RapidityCut = 1.0;
    
    // default
    float const cosTheta = 0.95;
    float const kDca[nCent][nPtBins] = {
        {0.0106, 0.0106, 0.0069, 0.0068, 0.0050, 0.0050}, // 60-80%
        {0.0140, 0.0100, 0.0075, 0.0072, 0.0060, 0.0050}, // 40-60%
        {0.0151, 0.0102, 0.0104, 0.0099, 0.0063, 0.0050}, // 20-40%
        {0.0145, 0.0113, 0.0094, 0.0089, 0.0069, 0.0050}, // 10-20%
        {0.0138, 0.0109, 0.0082, 0.0094, 0.0076, 0.0054}  // 0-10%
    };
    float const pDca[nCent][nPtBins] = {
        {0.0098, 0.0098, 0.0083, 0.0073, 0.0056, 0.0050}, // 60-80%
        {0.0145, 0.0128, 0.0072, 0.0079, 0.0060, 0.0051}, // 40-60%
        {0.0131, 0.0113, 0.0099, 0.0106, 0.0065, 0.0052}, // 20-40%
        {0.0141, 0.0100, 0.0074, 0.0077, 0.0066, 0.0052}, // 10-20%
        {0.0133, 0.0105, 0.0093, 0.0097, 0.0067, 0.0055}  // 0-10%
    };
    float const dcaV0ToPv[nCent][nPtBins] = {
        {0.0076, 0.0076, 0.0053, 0.0054, 0.0054, 0.0042}, // 60-80%
        {0.0072, 0.0057, 0.0058, 0.0049, 0.0049, 0.0047}, // 40-60%
        {0.0066, 0.0055, 0.0053, 0.0046, 0.0041, 0.0050}, // 20-40%
        {0.0063, 0.0047, 0.0045, 0.0046, 0.0042, 0.0044}, // 10-20%
        {0.0062, 0.0055, 0.0040, 0.0040, 0.0040, 0.0044}  // 0-10%
    };
    float const dcaDaughters[nCent][nPtBins] = {
        {0.0077, 0.0077, 0.0094, 0.0078, 0.0081, 0.0120}, // 60-80%
        {0.0080, 0.0083, 0.0092, 0.0081, 0.0094, 0.0106}, // 40-60%
        {0.0078, 0.0073, 0.0080, 0.0093, 0.0096, 0.0103}, // 20-40%
        {0.0076, 0.0078, 0.0092, 0.0072, 0.0086, 0.0085}, // 10-20%
        {0.0071, 0.0064, 0.0070, 0.0063, 0.0082, 0.0080}  // 0-10%
    };
    float const decayLength[nCent][nPtBins] = {
        {0.0175, 0.0175, 0.0187, 0.0178, 0.0184, 0.0187}, // 60-80%
        {0.0171, 0.0196, 0.0210, 0.0187, 0.0190, 0.0214}, // 40-60%
        {0.0178, 0.0206, 0.0221, 0.0209, 0.0219, 0.0240}, // 20-40%
        {0.0172, 0.0215, 0.0252, 0.0232, 0.0236, 0.0237}, // 10-20%
        {0.0100, 0.0199, 0.0227, 0.0232, 0.0236, 0.0255}  // 0-10%
    };
    
    // tight
    float const cosTheta1 = 0.95;
    float const kDca1[nCent][nPtBins] = {
        {0.0126, 0.0126, 0.0116, 0.0097, 0.0076, 0.0050}, // 60-80%
        {0.0176, 0.0100, 0.0123, 0.0092, 0.0068, 0.0056}, // 40-60%
        {0.0158, 0.0123, 0.0133, 0.0125, 0.0103, 0.0053}, // 20-40%
        {0.0172, 0.0165, 0.0119, 0.0125, 0.0091, 0.0057}, // 10-20%
        {0.0170, 0.0109, 0.0117, 0.0109, 0.0111, 0.0059}  // 0-10%
    };
    float const pDca1[nCent][nPtBins] = {
        {0.0130, 0.0130, 0.0130, 0.0095, 0.0097, 0.0086}, // 60-80%
        {0.0143, 0.0100, 0.0072, 0.0145, 0.0113, 0.0095}, // 40-60%
        {0.0150, 0.0113, 0.0063, 0.0149, 0.0107, 0.0105}, // 20-40%
        {0.0172, 0.0100, 0.0144, 0.0133, 0.0130, 0.0092}, // 10-20%
        {0.0139, 0.0148, 0.0093, 0.0133, 0.0080, 0.0088}  // 0-10%
    };
    float const dcaV0ToPv1[nCent][nPtBins] = {
        {0.0058, 0.0058, 0.0051, 0.0036, 0.0036, 0.0031}, // 60-80%
        {0.0056, 0.0045, 0.0050, 0.0037, 0.0029, 0.0026}, // 40-60%
        {0.0062, 0.0060, 0.0037, 0.0038, 0.0030, 0.0026}, // 20-40%
        {0.0056, 0.0049, 0.0039, 0.0036, 0.0033, 0.0028}, // 10-20%
        {0.0055, 0.0053, 0.0037, 0.0027, 0.0026, 0.0025}  // 0-10%
    };
    float const dcaDaughters1[nCent][nPtBins] = {
        {0.0076, 0.0076, 0.0078, 0.0095, 0.0076, 0.0093}, // 60-80%
        {0.0088, 0.0059, 0.0081, 0.0083, 0.0073, 0.0108}, // 40-60%
        {0.0061, 0.0043, 0.0070, 0.0089, 0.0066, 0.0070}, // 20-40%
        {0.0077, 0.0049, 0.0042, 0.0056, 0.0053, 0.0119}, // 10-20%
        {0.0077, 0.0044, 0.0047, 0.0073, 0.0060, 0.0061}  // 0-10%
    };
    float const decayLength1[nCent][nPtBins] = {
        {0.0203, 0.0203, 0.0206, 0.0228, 0.0161, 0.0216}, // 60-80%
        {0.0222, 0.0229, 0.0269, 0.0236, 0.0232, 0.0182}, // 40-60%
        {0.0240, 0.0242, 0.0268, 0.0319, 0.0176, 0.0338}, // 20-40%
        {0.0219, 0.0240, 0.0213, 0.0231, 0.0261, 0.0399}, // 10-20%
        {0.0100, 0.0230, 0.0268, 0.0292, 0.0249, 0.0225}  // 0-10%
    };
    
    // loose
    float const cosTheta2 = 0.95;
    float const kDca2[nCent][nPtBins] = {
        {0.0090, 0.0090, 0.0073, 0.0057, 0.0050, 0.0050}, // 60-80%
        {0.0115, 0.0114, 0.0067, 0.0050, 0.0050, 0.0050}, // 40-60%
        {0.0134, 0.0112, 0.0074, 0.0063, 0.0050, 0.0050}, // 20-40%
        {0.0135, 0.0120, 0.0070, 0.0060, 0.0050, 0.0050}, // 10-20%
        {0.0111, 0.0109, 0.0080, 0.0067, 0.0050, 0.0050}  // 0-10%
    };
    float const pDca2[nCent][nPtBins] = {
        {0.0098, 0.0098, 0.0080, 0.0058, 0.0050, 0.0050}, // 60-80%
        {0.0124, 0.0114, 0.0072, 0.0050, 0.0050, 0.0050}, // 40-60%
        {0.0111, 0.0113, 0.0089, 0.0062, 0.0050, 0.0050}, // 20-40%
        {0.0108, 0.0100, 0.0074, 0.0069, 0.0050, 0.0050}, // 10-20%
        {0.0116, 0.0105, 0.0093, 0.0072, 0.0050, 0.0050}  // 0-10%
    };
    float const dcaV0ToPv2[nCent][nPtBins] = {
        {0.0083, 0.0083, 0.0090, 0.0077, 0.0059, 0.0055}, // 60-80%
        {0.0069, 0.0083, 0.0087, 0.0077, 0.0085, 0.0072}, // 40-60%
        {0.0075, 0.0073, 0.0056, 0.0053, 0.0100, 0.0090}, // 20-40%
        {0.0071, 0.0063, 0.0055, 0.0052, 0.0065, 0.0055}, // 10-20%
        {0.0078, 0.0063, 0.0054, 0.0045, 0.0062, 0.0047}  // 0-10%
    };
    float const dcaDaughters2[nCent][nPtBins] = {
        {0.0097, 0.0097, 0.0123, 0.0092, 0.0093, 0.0092}, // 60-80%
        {0.0095, 0.0084, 0.0100, 0.0108, 0.0130, 0.0113}, // 40-60%
        {0.0081, 0.0090, 0.0078, 0.0094, 0.0137, 0.0121}, // 20-40%
        {0.0082, 0.0086, 0.0099, 0.0102, 0.0121, 0.0105}, // 10-20%
        {0.0098, 0.0087, 0.0088, 0.0078, 0.0100, 0.0101}  // 0-10%
    };
    float const decayLength2[nCent][nPtBins] = {
        {0.0154, 0.0154, 0.0163, 0.0147, 0.0126, 0.0140}, // 60-80%
        {0.0158, 0.0153, 0.0172, 0.0150, 0.0126, 0.0164}, // 40-60%
        {0.0100, 0.0177, 0.0177, 0.0194, 0.0131, 0.0150}, // 20-40%
        {0.0157, 0.0197, 0.0222, 0.0180, 0.0155, 0.0189}, // 10-20%
        {0.0148, 0.0179, 0.0215, 0.0198, 0.0167, 0.0206}  // 0-10%
    };

    TopologicalCuts const D0Cuts_50eff(RapidityCut,nPtBins,PtEdge,nCent,CentEdge,dcaV0ToPv1,decayLength1,dcaDaughters1,kDca1,pDca1);
    TopologicalCuts const D0Cuts_150eff(RapidityCut,nPtBins,PtEdge,nCent,CentEdge,dcaV0ToPv2,decayLength2,dcaDaughters2,kDca2,pDca2);
    TopologicalCuts const D0Cuts(RapidityCut,nPtBins,PtEdge,nCent,CentEdge,dcaV0ToPv,decayLength,dcaDaughters,kDca,pDca);
}
#endif
